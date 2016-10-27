/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2011, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */

#ifndef PCL_REGISTRATION_CORRESPONDENCE_ESTIMATION_LOOKUP_TABLE_H_
#define PCL_REGISTRATION_CORRESPONDENCE_ESTIMATION_LOOKUP_TABLE_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <pcl/pcl_base.h>
#include <pcl/pcl_macros.h>
#include <pcl/common/common.h>
#include <pcl/common/transforms.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/registration/correspondence_types.h>
#include <pcl/registration/correspondence_estimation.h>

namespace pcl
{
  namespace registration
  {
    /** \brief @b CorrespondenceLookupTableCell represents a match between
      * a CorrespondenceLookupTable cell and the closest point in the reference point cloud.
      *
      * Note:
      *   Use DistanceT=float  when initializing the lookup table using a k-d tree (the k-d tree returns floats for the distances, and as such, it should be used floats to save memory space).
      *   Use DistanceT=double when using the Meijster Squared Euclidean Distance Transform, because the Meijster algorithm uses squared integer distances (distance between two adjacent cells is 1, not cell_resolution)
      *   when initializing the lookup table and since it is undesirable the loss of precision to ensure correct functioning of the algorithm,
      *   it should be used doubles, since they can store integer numbers up to 2^53 while floats only store until 2^24 (maximum value stored in the mantissa).
      *
      * \author Carlos M. Costa
      * \ingroup registration
      */
    template <typename DistanceT = double>
    struct CorrespondenceLookupTableCell
    {
      /** \brief Index of the closest point in the reference point cloud. */
      int closest_point_index;

      /** \brief Squared distance to the closest point. */
      DistanceT distance_squared_to_closest_point;

      /** \brief Empty constructor.
        * Sets \ref closest_point_index to -1, \ref distance_to_closest_point to -FLT_MAX.
        */
      inline CorrespondenceLookupTableCell ()
      : closest_point_index (-1), distance_squared_to_closest_point (-std::numeric_limits<DistanceT>::max ()) {}

      /** \brief Constructor. */
      inline CorrespondenceLookupTableCell (int index, DistanceT distance)
      : closest_point_index (index), distance_squared_to_closest_point(distance) {}

      /** \brief Empty destructor. */
      virtual ~CorrespondenceLookupTableCell () {}
    };

    /** \brief @b LookupTable provides fast correspondence estimation
      * by pre-computing the closest point for each cell within a uniform volume grid.
      * It is recommended for 2D point cloud registration or 3D point cloud registration of
      * small volumes (it requires more memory than k-d trees but it is much faster).
      *
      * \author Carlos M. Costa
      * \ingroup registration
      */
    template <typename PointT, typename DistanceT = double>
    class CorrespondenceLookupTable
    {
      public:
        /** \brief Empty constructor. */
        CorrespondenceLookupTable ()
          : cell_resolution_ (0.01)
          , cell_resolution_inverse_ (100.0)
          , lookup_table_margin_ (1.0, 1.0, 1.0)
          , number_cells_x_ (0)
          , number_cells_y_ (0)
          , number_cells_z_ (0)
          , number_cells_xy_slice_ (0)
          , use_search_tree_when_query_point_is_outside_lookup_table_ (true)
          , compute_distance_from_query_point_to_closest_point_ (false)
          , initialize_lookup_table_using_euclidean_distance_transform_ (false)
          , number_of_queries_on_lookup_table_ (0)
          , number_of_queries_on_search_tree_ (0) {}

        /** \brief Empty destructor. */
        virtual ~CorrespondenceLookupTable () {}

        /** \brief Set the lookup table cell resolution.
          * \param[in] cell_resolution is the lookup table cell size.
          */
        inline void
        setCellResolution (float cell_resolution) { cell_resolution_ = cell_resolution; cell_resolution_inverse_ = 1.0 / cell_resolution_; }

        /** \brief Get the lookup table cell size. */
        inline float
        getCellResolution () { return (cell_resolution_); }

        /** \brief Set the lookup table margin (in x, y and z axis).
          * \param[in] lookup_table_margin is the extra space that will be filled with extra cells around the data bounding box
          * in order to have pre-computed correspondences available when registering point clouds that are not aligned.
          */
        inline void
        setLookupTableMargin (Eigen::Vector3f lookup_table_margin) { lookup_table_margin_ = lookup_table_margin; }

        /** \brief Gets the lookup table margin. */
        inline Eigen::Vector3f
        getLookupTableMargin () { return (lookup_table_margin_); }

        /** \brief Gets the lookup table minimum bounds. */
        inline Eigen::Vector3f
        getMinimumBounds () { return (minimum_bounds_); }

        /** \brief Gets the lookup table maximum bounds. */
        inline Eigen::Vector3f
        getMaximumBounds () { return (maximum_bounds_); }

        /** \brief Get the number of cells in the X axis. */
        inline size_t
        getNumberOfCellsXAxis () { return (number_cells_x_); }

        /** \brief Get the number of cells in the Y axis. */
        inline size_t
        getNumberOfCellsYAxis () { return (number_cells_y_); }

        /** \brief Get the number of cells in the Z axis. */
        inline size_t
        getNumberOfCellsZAxis () { return (number_cells_z_); }

        /** \brief Get a reference to the lookup table. */
        inline std::vector<pcl::registration::CorrespondenceLookupTableCell<DistanceT> >&
        getLookupTable () { return (lookup_table_); }

        /** \brief Get the search tree. */
        inline typename pcl::search::Search<PointT>::Ptr
        getSearchTree () { return (search_tree_); }

        /** \brief Set if the search tree shall be used as a fall back strategy when the query points are outside of the lookup table bounds.
          * \param[in] use_search_tree_when_query_point_is_outside_lookup_table True for using the search tree as a fall back strategy when the query points are outside the lookup table bounds.
          */
        inline void
        setUseSearchTreeWhenQueryPointIsOutsideLookupTable (bool use_search_tree_when_query_point_is_outside_lookup_table) { use_search_tree_when_query_point_is_outside_lookup_table_ = use_search_tree_when_query_point_is_outside_lookup_table; }

        /** \brief Check if the search tree will be used as a fall back strategy when the query points are outside of the lookup table bounds. */
        inline bool
        getUseSearchTreeWhenQueryPointIsOutsideLookupTable () { return (use_search_tree_when_query_point_is_outside_lookup_table_); }

        /** \brief Set if the distance between query and closest point should be computed.
          * \param[in] compute_distance_from_query_point_to_closest_point True for computing the distance between query point and the closest point. False for using the distance between the centroids of the cells associated with the query and closest point.
          */
        inline void
        setComputeDistanceFromQueryPointToClosestPoint (bool compute_distance_from_query_point_to_closest_point) { compute_distance_from_query_point_to_closest_point_ = compute_distance_from_query_point_to_closest_point; }

        /** \brief Check if the distance between query and closest point will be computed. */
        inline bool
        getComputeDistanceFromQueryPointToClosestPoint () { return (compute_distance_from_query_point_to_closest_point_); }

        /** \brief Set the lookup table initialization algorithm.
          * \param[in] initialize_lookup_table_using_euclidean_distance_transform True for using the Euclidean Distance Transform (much faster). False for using a k-d tree (more accurate).
          */
        inline void
        setInitializeLookupTableUsingEuclideanDistanceTransform (bool initialize_lookup_table_using_euclidean_distance_transform) { initialize_lookup_table_using_euclidean_distance_transform_ = initialize_lookup_table_using_euclidean_distance_transform; }

        /** \brief Get the lookup table initialization method (true if using the Euclidean Distance Transform or false if relying on a k-d tree). */
        inline bool
        getInitializeLookupTableUsingEuclideanDistanceTransform () { return (initialize_lookup_table_using_euclidean_distance_transform_); }

        /** \brief Gets the number of queries performed on the lookup table. */
        inline size_t
        getNumberOfQueriesOnLookupTable () { return (number_of_queries_on_lookup_table_); }

        /** \brief Resets the number of queries performed on the lookup table. */
        inline void
        resetNumberOfQueriesOnLookupTable () { number_of_queries_on_lookup_table_ = 0; }

        /** \brief Gets the number of queries performed on the lookup table. */
        inline size_t
        getNumberOfQueriesOnSearchTree () { return (number_of_queries_on_search_tree_); }

        /** \brief Resets the number of queries performed on the lookup table. */
        inline void
        resetNumberOfQueriesOnSearchTree () { number_of_queries_on_search_tree_ = 0; }

        /**
          * \brief Computes the lookup table minimum and maximum bounds with the additional margin that was set previously.
          * @param[in] pointcloud Point cloud with the data.
          * @return True if the computed bounds are valid.
          */
        virtual bool
        computeLookupTableBounds (const pcl::PointCloud<PointT>& pointcloud);

        /**
          * \brief Initialize the lookup table using the provided tree.
          * @param[in] tree Search tree with the point data.
          * @return True if the lookup table was initialized successfully.
          */
        virtual bool
        initLookupTable (typename pcl::search::Search<PointT>::ConstPtr tree);

        /**
          * \brief Initialize the lookup table using the Euclidean Distance Transform algorithm.
          *
          * \note Implementation based on the algorithm introduced in:
          * A. Meijster, J. B. T. M. Roerdink and W. H. Hesselink, A general algorithm for computing distance transforms in linear time.
          * In: Mathematical Morphology and its Applications to Image and Signal Processing, Kluwer Acad. Publ., 2000, pp. 331-340.
          *
          * @param[in] pointcloud Point cloud with the data.
          * @param[out] discarded_indices The distance transform uses the voxel centroids when computing distances.
          *   This means that if several reference points fall within the same voxel, only the point closest to the voxel centroid will be kept.
          *   If discarded_indices != null, then it will be retrieved which point indices were discarded (because there was another point within the same voxel and closest to the centroid).
          * @return True if the lookup table was initialized successfully.
          */
        virtual bool
        initLookupTableUsingEuclideanDistanceTransform (typename pcl::PointCloud<PointT>::ConstPtr pointcloud, pcl::IndicesPtr discarded_indices);

        /**
         * Gets the correspondence cell index associated with the query_point.
         * @param[in] query_point Coordinates of the query point.
         * @param[out] correspondence_index The computed index of the lookup table cell.
         * @param[out] correspondence_index_components The computed index [x,y,z] components.
         * @return True if the index is valid. False if the query point was outside of the lookup table bounds.
         */
        bool
        computeCorrespondenceCellIndex (const PointT& query_point, size_t& correspondence_index, Eigen::Vector3i& correspondence_index_components);

        /**
         * Computes the correspondence cell centroid from the lookup table cell index components (retrieved with @see computeCorrespondenceCellIndex).
         * @return The centroid of the correspondence cell.
         */
        inline Eigen::Vector3f
        computeCorrespondenceCellCentroid (const Eigen::Vector3i& correspondence_index_components)
        {
          return (Eigen::Vector3f (minimum_bounds_ (0) + (correspondence_index_components (0) * cell_resolution_),
                                   minimum_bounds_ (1) + (correspondence_index_components (1) * cell_resolution_),
                                   minimum_bounds_ (2) + (correspondence_index_components (2) * cell_resolution_)));
        }

        /**
         * Computes the squared distance between two points.
         * @param a First point.
         * @param b Second point.
         * @return The Euclidean squared distance between the two points.
         */
        inline float
        computeSquaredDistance (const Eigen::Vector3f& a, const Eigen::Vector3f& b)
        {
          return ((a - b).squaredNorm ());
        }

        /**
         * Computes the squared distance between two points.
         * @param a First point.
         * @param b Second point.
         * @return The Euclidean squared distance between the two points.
         */
        inline float
        computeSquaredDistance (const PointT& a, const PointT& b)
        {
          float dx = a.x - b.x;
          float dy = a.y - b.y;
          float dz = a.z - b.z;
          return (dx * dx + dy * dy + dz * dz);
        }

        /**
         * Computes the squared distance between the query point and the cell centroid.
         * @param query_point Query point.
         * @param centroid Cell centroid.
         * @return The Euclidean squared distance between the two points.
         */
        inline float
        computeSquaredDistanceToCentroid (const PointT& query_point, Eigen::Vector3f centroid)
        {
          float dx = query_point.x - centroid (0);
          float dy = query_point.y - centroid (1);
          float dz = query_point.z - centroid (2);
          return (dx * dx + dy * dy + dz * dz);
        }

        /**
         * Retrieves the point associated to a given cell by using the index of the point cloud that was used to initialize the lookup table.
         * @param[in] index Index of the point to retrieve.
         * @param[out] point Point retrieved.
         * @return True if the index was valid and the point was retrieved.
         */
        inline bool
        getCellPoint (size_t index, PointT& point)
        {
          if (pointcloud_ && index < pointcloud_->size ())
          {
            point = (*pointcloud_)[index];
            return (true);
          }
          return (false);
        }

        /**
         * Gets the pre-computed correspondence associated with the provided query point.
         * @param[in] query_point Coordinates of the query point.
         * @param[in] maximum_correspondence_distance_squared The maximum distance squared that a valid correspondence can have.
         * @param[out] correspondence The correspondence found for the provided query point. If the query point is outside the pre-computed cells,
         * the search tree will be used for finding the correspondence.
         * @return True if it was found a valid correspondence.
         */
        bool
        getCorrespondence (const PointT& query_point, double maximum_correspondence_distance_squared, pcl::registration::CorrespondenceLookupTableCell<DistanceT>& correspondance);

        /**
         * Gets the pre-computed correspondence using a valid index computed with @see computeCorrespondenceCellIndex.
         * @param correspondence_index The index of the lookup table cell.
         * @return Reference to the correspondence on the provided index.
         */
        inline pcl::registration::CorrespondenceLookupTableCell<DistanceT>&
        getCorrespondence (size_t correspondence_index)
        {
          return (lookup_table_[correspondence_index]);
        }

        /**
         * Gets the pre-computed correspondence associated with the provided coordinates.
         * @param x x-axis cell coordinate
         * @param y y-axis cell coordinate
         * @param z z-axis cell coordinate
         * @return Reference to the correspondence.
         */
        inline pcl::registration::CorrespondenceLookupTableCell<DistanceT>&
        getCorrespondence (size_t x, size_t y, size_t z)
        {
          return (lookup_table_[x + y * number_cells_x_ + z * number_cells_xy_slice_]);
        }

      protected:
        /**
         * \brief Adds points to the lookup table. If several points fall within the same voxel, only the one closest to the voxel centroid will be kept.
         * @param[in] pointcloud Point cloud with the data.
         * @param[out] discarded_indices If discarded_indices != null, then it will be retrieved which point indices were discarded (because there was another point within the same voxel and closest to the centroid).
         */
        void
        addPointsToLookupTable (const pcl::PointCloud<PointT>& pointcloud, pcl::IndicesPtr discarded_indices);

        /** \brief Computes the first phase of the Meijster Euclidean distance transform algorithm on the data x axis */
        void
        computeFirstPhaseMeijsterAxisX ();

        /** \brief Computes the second phase of the Meijster Euclidean distance transform algorithm on the data y axis */
        void
        computeSecondPhaseMeijsterAxisY ();

        /** \brief Computes the second phase of the Meijster Euclidean distance transform algorithm on the data z axis */
        void
        computeSecondPhaseMeijsterAxisZ ();

        /** \brief Updates the distance within each lookup table cell to take in consideration the cell resolution.
         *
         * The Meijster algorithm assumes that the distance between two adjacent cells is one unit.
         * This method scales the distance results of the Meijster algorithm in order to have distance_between_adjacent_cells = cell_resolution. */
        void
        updateCellDistancesWithCellResolution ();

        /**
         * Computes Meijster squared Euclidean distance metric.
         * @param x Current x cell index.
         * @param i Index of the closest interest cell (in a given axis).
         * @param gi_squared Distance squared to the closest interest cell stored in the i cell.
         * @return Squared Euclidean distance.
         */
        inline float
        computeMeijsterSquaredEDT (size_t x, size_t i, DistanceT gi_squared)
        {
          int xi = (x - i);
          return (xi * xi + gi_squared);
        }

        /**
         * Computes the Meijster separator index (Sep) corresponding to the iu parabola intersection.
         * @param i Index of the closest interest cell (in a given axis).
         * @param u Upper bound of i index.
         * @param gi_squared Distance to the closest interest cell stored in the i cell.
         * @param gu_squared Distance to the closest interest cell stored in the u cell.
         * @return Meijster Sep value.
         */
        inline size_t
        computeMeijsterParabolaIntersectionIndex (size_t i, size_t u, DistanceT gi_squared, DistanceT gu_squared)
        {
          // (x - i)^2 + (gi)^2 <= (x - u)^2 + (gu)^2
          // x <= (u^2 - i^2 + (gu)^2 - (gi)^2) / (2(u-i))
          return ((u * u - i * i + gu_squared - gi_squared) / (2 * (u - i)));
        }

        /** \brief 3 Dimensional array containing the pre-computed correspondences. */
        std::vector<pcl::registration::CorrespondenceLookupTableCell<DistanceT> > lookup_table_;

        /** \brief Search tree associated with the lookup table. */
        typename pcl::search::Search<PointT>::ConstPtr search_tree_;

        /** \brief Point cloud associated with the lookup table. */
        typename pcl::PointCloud<PointT>::ConstPtr pointcloud_;

        /** \brief Resolution of the lookup table. */
        float cell_resolution_;

        /** \brief Inverse of the resolution of the lookup table. */
        float cell_resolution_inverse_;

        /** \brief Margin added to the data bounding box in order to create extra cells in the lookup table surrounding the point cloud data. */
        Eigen::Vector3f lookup_table_margin_;

        /** \brief Lookup table minimum bounds. */
        Eigen::Vector3f minimum_bounds_;

        /** \brief Lookup table maximum bounds. */
        Eigen::Vector3f maximum_bounds_;

        /** \brief Number of cells in the x dimension. */
        size_t number_cells_x_;

        /** \brief Number of cells in the y dimension. */
        size_t number_cells_y_;

        /** \brief Number of cells in the z dimension. */
        size_t number_cells_z_;

        /** \brief Precomputed number of cells for each XY plane / slice */
        size_t number_cells_xy_slice_;

        /** \brief If true, the search tree will be used as a fall back strategy when the query points are outside of the lookup table bounds. */
        bool use_search_tree_when_query_point_is_outside_lookup_table_;

        /** \brief The lookup table stores the distance between the cell centroids.
         * If it is required accurate distance estimation, this flag should be set to true for computing the distance between each query point and its correspondence. */
        bool compute_distance_from_query_point_to_closest_point_;

        /** \brief If true it will be used the Euclidean Distance Transform for initializing the lookup table. Otherwise it will be used a k-d tree. */
        bool initialize_lookup_table_using_euclidean_distance_transform_;

        /** \brief Number of queries performed on the lookup table. */
        size_t number_of_queries_on_lookup_table_;

        /** \brief Number of queries performed on the search tree (because they were outside the lookup table bounds). */
        size_t number_of_queries_on_search_tree_;
    };

    /** \brief @b CorrespondenceEstimationLookupTable provides fast correspondence estimation
      * by using a lookup table with pre-computed correspondences.
      *
      * \author Carlos M. Costa
      * \ingroup registration
      */
    template <typename PointSource, typename PointTarget, typename Scalar = float, typename DistanceT = double>
    class CorrespondenceEstimationLookupTable: public CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>
    {
      public:
        typedef boost::shared_ptr<CorrespondenceEstimationLookupTable<PointSource, PointTarget, Scalar, DistanceT> > Ptr;
        typedef boost::shared_ptr<const CorrespondenceEstimationLookupTable<PointSource, PointTarget, Scalar, DistanceT> > ConstPtr;

        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::corr_name_;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::force_no_recompute_;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::force_no_recompute_reciprocal_;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::indices_;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::initCompute;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::initComputeReciprocal;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::input_;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::target_;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::tree_;
        using CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::tree_reciprocal_;
        using PCLBase<PointSource>::deinitCompute;

        typedef pcl::search::KdTree<PointTarget> KdTree;
        typedef typename pcl::search::KdTree<PointTarget>::Ptr KdTreePtr;

        typedef pcl::search::KdTree<PointSource> KdTreeReciprocal;
        typedef typename KdTree::Ptr KdTreeReciprocalPtr;

        typedef pcl::PointCloud<PointSource> PointCloudSource;
        typedef typename PointCloudSource::Ptr PointCloudSourcePtr;
        typedef typename PointCloudSource::ConstPtr PointCloudSourceConstPtr;

        typedef pcl::PointCloud<PointTarget> PointCloudTarget;
        typedef typename PointCloudTarget::Ptr PointCloudTargetPtr;
        typedef typename PointCloudTarget::ConstPtr PointCloudTargetConstPtr;

        typedef typename KdTree::PointRepresentationConstPtr PointRepresentationConstPtr;

        /** \brief Empty constructor. */
        CorrespondenceEstimationLookupTable ()
        {
          corr_name_  = "CorrespondenceEstimationLookupTable";
        }

        /** \brief Empty destructor. */
        virtual ~CorrespondenceEstimationLookupTable () {}

        /** \brief Get the source lookup table. */
        inline CorrespondenceLookupTable<PointSource, DistanceT>&
        getSourceCorrespondencesLookupTable () { return (source_correspondences_lookup_table_); }

        /** \brief Get the target lookup table. */
        inline CorrespondenceLookupTable<PointTarget, DistanceT>&
        getTargetCorrespondencesLookupTable () { return (target_correspondences_lookup_table_); }

        /** \brief Internal computation initialization. */
        bool
        initComputeReciprocal ();

        /** \brief Internal computation initialization for reciprocal correspondences. */
        bool
        initCompute ();

        /** \brief Provide a pointer to the search object used for finding correspondences in
          * the source cloud (used for reciprocal correspondence finding).
          * \param[in] tree a pointer to the spatial search object.
          * \param[in] force_no_recompute If set to true, this tree and the associated lookup table
          * will never be recomputed, regardless of calls to setInputSource.
          */
        virtual void
        setSearchMethodSource (const KdTreeReciprocalPtr &tree,
                               bool force_no_recompute = false)
        {
          CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::setSearchMethodSource (tree, force_no_recompute);
          //if (!force_no_recompute)
          //  source_correspondences_lookup_table_.initLookupTable(tree);
        }

        /** \brief Provide a pointer to the search object used for finding correspondences in the target cloud.
          * \param[in] tree a pointer to the spatial search object.
          * \param[in] force_no_recompute If set to true, this tree and the associated lookup table
          * will never be recomputed, regardless of calls to setInputTarget.
          */
        virtual void
        setSearchMethodTarget (const KdTreePtr &tree,
                               bool force_no_recompute = false)
        {
          CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::setSearchMethodTarget (tree, force_no_recompute);
          if (!force_no_recompute)
            target_correspondences_lookup_table_.initLookupTable(tree);
        }

        /** \brief Determine the correspondences between input and target cloud.
          * \param[out] correspondences The found correspondences (index of query point, index of target point, distance).
          * \param[in] max_distance Maximum allowed distance between correspondences.
          */
        virtual void
        determineCorrespondences (pcl::Correspondences &correspondences,
                                  double max_distance = std::numeric_limits<double>::max ());

        /** \brief Determine the reciprocal correspondences between input and target cloud.
          * A correspondence is considered reciprocal if the Src_i -> Tgt_i correspondence is the same as Tgt_i -> Src_i.
          *
          * \param[out] correspondences The found correspondences (index of query and target point, distance).
          * \param[in] max_distance Maximum allowed distance between correspondences.
          */
        virtual void
        determineReciprocalCorrespondences (pcl::Correspondences &correspondences,
                                            double max_distance = std::numeric_limits<double>::max ());

        /** \brief Clone and cast to CorrespondenceEstimationBase. */
        virtual boost::shared_ptr< CorrespondenceEstimationBase<PointSource, PointTarget, Scalar> >
        clone () const
        {
          Ptr copy (new CorrespondenceEstimationLookupTable<PointSource, PointTarget, Scalar, DistanceT> (*this));
          return (copy);
        }

      protected:
        /** \brief LookupTable containing the pre-computed source correspondences */
        CorrespondenceLookupTable<PointSource, DistanceT> source_correspondences_lookup_table_;

        /** \brief LookupTable containing the pre-computed target correspondences */
        CorrespondenceLookupTable<PointTarget, DistanceT> target_correspondences_lookup_table_;
     };
  }
}

#include <pcl/registration/impl/correspondence_estimation_lookup_table.hpp>

#endif /* PCL_REGISTRATION_CORRESPONDENCE_ESTIMATION_LOOKUP_TABLE_H_ */
