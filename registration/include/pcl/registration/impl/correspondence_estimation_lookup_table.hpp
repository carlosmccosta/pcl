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

#ifndef PCL_REGISTRATION_IMPL_CORRESPONDENCE_ESTIMATION_LOOKUP_TABLE_H_
#define PCL_REGISTRATION_IMPL_CORRESPONDENCE_ESTIMATION_LOOKUP_TABLE_H_

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
bool
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::computeLookupTableBounds (const pcl::PointCloud<PointT>& pointcloud)
{
  if (pointcloud.size () < 2)
    return (false);

  Eigen::Vector4f min_pt;
  Eigen::Vector4f max_pt;
  pcl::getMinMax3D (pointcloud, min_pt, max_pt);
  float half_cell_resolution = cell_resolution_ * 0.5;
  for (size_t i = 0; i < 3; ++i)
  {
    minimum_bounds_ (i) = min_pt (i) - (lookup_table_margin_ (i) + half_cell_resolution);
    maximum_bounds_ (i) = max_pt (i) + (lookup_table_margin_ (i) + half_cell_resolution);
  }

  number_cells_x_ = std::max ((size_t)(std::ceil ((maximum_bounds_ (0) - minimum_bounds_ (0)) * cell_resolution_inverse_)), (size_t)1);
  number_cells_y_ = std::max ((size_t)(std::ceil ((maximum_bounds_ (1) - minimum_bounds_ (1)) * cell_resolution_inverse_)), (size_t)1);
  number_cells_z_ = std::max ((size_t)(std::ceil ((maximum_bounds_ (2) - minimum_bounds_ (2)) * cell_resolution_inverse_)), (size_t)1);
  number_cells_xy_slice_ = number_cells_x_ * number_cells_y_; // precompute the number of cells within each xy slice.

  return (true);
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
bool
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::initLookupTable (typename pcl::search::Search<PointT>::ConstPtr tree)
{
  if (initialize_lookup_table_using_euclidean_distance_transform_)
  {
    if (initLookupTableUsingEuclideanDistanceTransform (tree->getInputCloud (), pcl::IndicesPtr ()))
    {
      search_tree_ = tree;
      return true;
    }
    return false;
  }

  if (!computeLookupTableBounds (*(tree->getInputCloud ())))
    return (false);

  search_tree_ = tree;
  pointcloud_ = tree->getInputCloud ();
  size_t number_cells = number_cells_z_ * number_cells_y_ * number_cells_x_;
  lookup_table_.clear ();
  lookup_table_.resize (number_cells);

  std::vector<int> index (1);
  std::vector<float> distance (1);
  int correspondence_number = 0;
  float half_cell_resolution = cell_resolution_ * 0.5;
  PointT query_point;
  query_point.z = minimum_bounds_(2) + half_cell_resolution;

  for (size_t z = 0; z < number_cells_z_; ++z)
  {
    query_point.y = minimum_bounds_(1) + half_cell_resolution;
    for (size_t y = 0; y < number_cells_y_; ++y)
    {
      query_point.x = minimum_bounds_(0) + half_cell_resolution;
      for (size_t x = 0; x < number_cells_x_; ++x)
      {
        search_tree_->nearestKSearch (query_point, 1, index, distance);
        pcl::registration::CorrespondenceLookupTableCell<DistanceT>& correspondence = lookup_table_[correspondence_number++];
        correspondence.closest_point_index = index[0];
        correspondence.distance_squared_to_closest_point = distance[0];
        query_point.x += cell_resolution_;
      }
      query_point.y += cell_resolution_;
    }
    query_point.z += cell_resolution_;
  }

  PCL_DEBUG ("[pcl::registration::LookupTable::initLookupTable] Initialized correspondence estimation LookupTable with %f resolution containing [x:%u|y:%u|z:%u]=%u cells from a point cloud with %u points\n",
             cell_resolution_, number_cells_x_, number_cells_y_, number_cells_z_, number_cells, tree->getInputCloud ()->size ());
  return (true);
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
bool
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::initLookupTableUsingEuclideanDistanceTransform (typename pcl::PointCloud<PointT>::ConstPtr pointcloud, pcl::IndicesPtr discarded_indices)
{
  if (!computeLookupTableBounds (*pointcloud))
    return (false);

  pointcloud_ = pointcloud;
  size_t number_cells = number_cells_z_ * number_cells_y_ * number_cells_x_;
  lookup_table_.clear ();
  lookup_table_.resize (number_cells);

  addPointsToLookupTable (*pointcloud_, discarded_indices);
  computeFirstPhaseMeijsterAxisX ();
  computeSecondPhaseMeijsterAxisY ();
  computeSecondPhaseMeijsterAxisZ ();
  updateCellDistancesWithCellResolution ();

  PCL_ERROR ("[pcl::registration::LookupTable::initLookupTableUsingEuclideanDistanceTransform] Initialized correspondence estimation LookupTable with %f resolution containing [x:%i|y:%i|z:%i]=%i cells from a point cloud with %u points\n",
             cell_resolution_, number_cells_x_, number_cells_y_, number_cells_z_, number_cells, pointcloud_->size ());
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
void
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::addPointsToLookupTable (const pcl::PointCloud<PointT>& pointcloud, pcl::IndicesPtr discarded_indices)
{
  size_t correspondence_index = 0;
  Eigen::Vector3i correspondence_index_components;
  for (size_t i = 0; i < pointcloud.size (); ++i)
  {
    if (computeCorrespondenceCellIndex (pointcloud[i], correspondence_index, correspondence_index_components))
    {
      pcl::registration::CorrespondenceLookupTableCell<DistanceT>& cell = getCorrespondence (correspondence_index);
      if (cell.closest_point_index < 0)
      {
        // First point to fall within the cell.
        cell.closest_point_index = i;
        cell.distance_squared_to_closest_point = 0;
      }
      else
      {
        Eigen::Vector3f cell_centroid = computeCorrespondenceCellCentroid (correspondence_index_components);
        float distance_of_existing_cell_point_to_centroid = computeSquaredDistanceToCentroid (pointcloud[cell.closest_point_index], cell_centroid);
        float distance_of_new_cell_point_to_centroid = computeSquaredDistanceToCentroid (pointcloud[i], cell_centroid);

        // Cell already has a point.
        // If the new point is closer to the centroid it will replace the current point.
        if (distance_of_new_cell_point_to_centroid < distance_of_existing_cell_point_to_centroid)
        {
          if (discarded_indices)
            discarded_indices->push_back (cell.closest_point_index);
          cell.closest_point_index = i;
        }
      }
    }
    else if (discarded_indices)
    {
      // Point is outside the lookup table bounds.
      discarded_indices->push_back (i);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
void
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::computeFirstPhaseMeijsterAxisX ()
{
  int number_cells_y_int = (int)number_cells_y_;
  for (size_t z = 0; z < number_cells_z_; ++z)
  {
    #pragma omp parallel for
    for (int y = 0; y < number_cells_y_int; ++y) // OpenMP supports size_t index only after version 3.0.
    {
      // Forward scan computes the integer distances to the closest interest cells (after they are found).
      size_t cell_index = 1 + (size_t)y * number_cells_x_ + z * number_cells_xy_slice_;
      for (size_t x = 1; x < number_cells_x_; ++x)
      {
        pcl::registration::CorrespondenceLookupTableCell<DistanceT>& previous_cell = getCorrespondence (cell_index - 1);
        pcl::registration::CorrespondenceLookupTableCell<DistanceT>& current_cell = getCorrespondence (cell_index++);
        if (previous_cell.distance_squared_to_closest_point >= 0)
        {
          current_cell.closest_point_index = previous_cell.closest_point_index;
          current_cell.distance_squared_to_closest_point = previous_cell.distance_squared_to_closest_point + 1;
        }
      }

      // Backwards scan propagates the integer distances of the found interest cells.
      --cell_index;
      for (int x = (int)number_cells_x_ - 2; x >= 0; --x)
      {
        pcl::registration::CorrespondenceLookupTableCell<DistanceT>& current_cell = getCorrespondence (cell_index);
        pcl::registration::CorrespondenceLookupTableCell<DistanceT>& next_cell = getCorrespondence (cell_index-- + 1);
        if (current_cell.closest_point_index < 0 || next_cell.distance_squared_to_closest_point < current_cell.distance_squared_to_closest_point)
        {
          current_cell.closest_point_index = next_cell.closest_point_index;
          current_cell.distance_squared_to_closest_point = next_cell.distance_squared_to_closest_point + 1;
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
void
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::computeSecondPhaseMeijsterAxisY ()
{
  int number_cells_x_int = (int)number_cells_x_;
  for (size_t z = 0; z < number_cells_z_; ++z)
  {
    #pragma omp parallel for
    for (int x = 0; x < number_cells_x_int; ++x) // OpenMP supports size_t index only after version 3.0.
    {
      // Forward scan computes the lower envelop segments.
      size_t cell_base_index = x + z * number_cells_xy_slice_;
      std::vector <size_t> y_index_where_closest_interest_cell_was_found (number_cells_y_, 0); // s
      std::vector <size_t> y_index_of_intersection_with_previous_parabola (number_cells_y_, 0); // t
      size_t last_segment_index = 0; // q
      DistanceT last_segment_index_distance_squared_to_closest_point = std::pow (getCorrespondence (cell_base_index).distance_squared_to_closest_point, 2);

      for (size_t y = 1; y < number_cells_y_; ++y)
      {
        DistanceT current_cell_distance_squared_to_closest_point = std::pow (getCorrespondence (cell_base_index + y * number_cells_x_).distance_squared_to_closest_point, 2);
        while (last_segment_index >= 0 &&
            computeMeijsterSquaredEDT (y_index_of_intersection_with_previous_parabola[last_segment_index],
                                       y_index_where_closest_interest_cell_was_found[last_segment_index],
                                       last_segment_index_distance_squared_to_closest_point) >
            computeMeijsterSquaredEDT (y_index_of_intersection_with_previous_parabola[last_segment_index],
                                       y,
                                       current_cell_distance_squared_to_closest_point))
        {
          // Fig. 3 (c) or Fig. 3 (b) cases (if q drops bellow 0).
          // Fy either intersects or is below the current lower envelop.
          // The current segments will either shrink or disappear.
          --last_segment_index;
          last_segment_index_distance_squared_to_closest_point = std::pow (getCorrespondence (cell_base_index + y_index_where_closest_interest_cell_was_found[last_segment_index] * number_cells_x_).distance_squared_to_closest_point, 2);
        }

        if (last_segment_index < 0)
        {
          // Fig. 3 (b) case.
          // Fy is below the current lower envelop. All old segments are replaced by a single new segment.
          last_segment_index = 0;
          last_segment_index_distance_squared_to_closest_point = std::pow (getCorrespondence (cell_base_index + y_index_where_closest_interest_cell_was_found[last_segment_index] * number_cells_x_).distance_squared_to_closest_point, 2);
          y_index_where_closest_interest_cell_was_found[0] = y;
        }
        else
        {
          // Fig. 3 (a) case if w >= number_cells_y_.
          // Fy is above the current lower envelop. It is discarded since it does not reduce the SEDT.
          float sep_plus_one = 1 + computeMeijsterParabolaIntersectionIndex (y_index_where_closest_interest_cell_was_found[last_segment_index],
                                                                             y,
                                                                             last_segment_index_distance_squared_to_closest_point,
                                                                             current_cell_distance_squared_to_closest_point);
          if (sep_plus_one < number_cells_y_)
          {
            // Fig. 3 (c) case.
            // Fy intersects the current lower envelop. A new segment is added.
            ++last_segment_index;
            last_segment_index_distance_squared_to_closest_point = std::pow (getCorrespondence (cell_base_index + y_index_where_closest_interest_cell_was_found[last_segment_index] * number_cells_x_).distance_squared_to_closest_point, 2);
            y_index_where_closest_interest_cell_was_found[last_segment_index] = y;
            y_index_of_intersection_with_previous_parabola[last_segment_index] = sep_plus_one;
          }
        }
      }

      // Backwards scan computes the squared distances.
      for (int y = (int)number_cells_y_ - 1; y >= 0; --y)
      {
        pcl::registration::CorrespondenceLookupTableCell<DistanceT>& current_cell = getCorrespondence (cell_base_index + y * number_cells_x_);
        current_cell.distance_squared_to_closest_point = computeMeijsterSquaredEDT (y, y_index_where_closest_interest_cell_was_found[last_segment_index], last_segment_index_distance_squared_to_closest_point);
        if (y == (int)y_index_of_intersection_with_previous_parabola[last_segment_index])
        {
          --last_segment_index;
          last_segment_index_distance_squared_to_closest_point = std::pow (getCorrespondence (cell_base_index + y_index_where_closest_interest_cell_was_found[last_segment_index] * number_cells_x_).distance_squared_to_closest_point, 2);
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
void
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::computeSecondPhaseMeijsterAxisZ ()
{
  // todo: finish second phase of Meijster in the Z axis
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
void
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::updateCellDistancesWithCellResolution ()
{
  int number_cells_y_int = (int)number_cells_y_;
  float cell_resolution_squared = cell_resolution_ * cell_resolution_;
  for (size_t z = 0; z < number_cells_z_; ++z)
  {
    #pragma omp parallel for
    for (int y = 0; y < number_cells_y_int; ++y) // OpenMP supports size_t index only after version 3.0.
    {
      size_t cell_index = (size_t)y * number_cells_x_ + z * number_cells_xy_slice_;
      for (size_t x = 0; x < number_cells_x_; ++x)
      {
        getCorrespondence (cell_index++).distance_squared_to_closest_point *= cell_resolution_squared;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
bool
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::computeCorrespondenceCellIndex (const PointT& query_point, size_t& correspondence_index, Eigen::Vector3i& correspondence_index_components)
{
  correspondence_index_components (0) = (int)((query_point.x - minimum_bounds_ (0)) * cell_resolution_inverse_);
  correspondence_index_components (1) = (int)((query_point.y - minimum_bounds_ (1)) * cell_resolution_inverse_);
  correspondence_index_components (2) = (int)((query_point.z - minimum_bounds_ (2)) * cell_resolution_inverse_);

  size_t index = correspondence_index_components (0) + correspondence_index_components (1) * number_cells_x_ + correspondence_index_components (2) * number_cells_xy_slice_;
  if (correspondence_index_components (0) >= 0 && (size_t)correspondence_index_components (0) < number_cells_x_ &&
      correspondence_index_components (1) >= 0 && (size_t)correspondence_index_components (1) < number_cells_y_ &&
      correspondence_index_components (2) >= 0 && (size_t)correspondence_index_components (2) < number_cells_z_ &&
      index < lookup_table_.size())
  {
    correspondence_index = index;
    return (true);
  }
  return (false);
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename DistanceT>
bool
pcl::registration::CorrespondenceLookupTable<PointT, DistanceT>::getCorrespondence (const PointT& query_point, double maximum_correspondence_distance_squared, pcl::registration::CorrespondenceLookupTableCell<DistanceT>& correspondance)
{
  size_t correspondence_index;
  Eigen::Vector3i correspondence_index_components;
  if (computeCorrespondenceCellIndex (query_point, correspondence_index, correspondence_index_components))
  {
    correspondance = lookup_table_[correspondence_index];
    ++number_of_queries_on_lookup_table_;

    if (compute_distance_from_query_point_to_closest_point_ && pointcloud_)
    {
      PointT closest_point;
      if (getCellPoint (correspondance.closest_point_index, closest_point))
        correspondance.distance_squared_to_closest_point = computeSquaredDistance (query_point, closest_point);
    }
  }
  else if (use_search_tree_when_query_point_is_outside_lookup_table_ && search_tree_)
  {
    std::vector<int> index (1);
    std::vector<float> distance (1);
    search_tree_->nearestKSearch (query_point, 1, index, distance);
    correspondance.closest_point_index = index[0];
    correspondance.distance_squared_to_closest_point = distance[0];
    ++number_of_queries_on_search_tree_;
  }
  else
    return (false);

  if (correspondance.distance_squared_to_closest_point <= maximum_correspondence_distance_squared)
    return (true);
  else
    return (false);
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar, typename DistanceT> bool
pcl::registration::CorrespondenceEstimationLookupTable<PointSource, PointTarget, Scalar, DistanceT>::initComputeReciprocal ()
{
  if (CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::initComputeReciprocal ())
  {
    if (!force_no_recompute_reciprocal_)
      return (source_correspondences_lookup_table_.initLookupTable(tree_reciprocal_));
    return (true);
  }
  return (false);
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar, typename DistanceT> bool
pcl::registration::CorrespondenceEstimationLookupTable<PointSource, PointTarget, Scalar, DistanceT>::initCompute ()
{
  if (CorrespondenceEstimationBase<PointSource, PointTarget, Scalar>::initCompute ())
  {
    if (!force_no_recompute_)
      return (target_correspondences_lookup_table_.initLookupTable(tree_));
    return (true);
  }
  return (false);
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar, typename DistanceT> void
pcl::registration::CorrespondenceEstimationLookupTable<PointSource, PointTarget, Scalar, DistanceT>::determineCorrespondences (
    pcl::Correspondences &correspondences, double max_distance)
{
  if (!initCompute ())
    return;

  double max_distance_squared = max_distance * max_distance;
  correspondences.resize (indices_->size ());
  pcl::registration::CorrespondenceLookupTableCell<DistanceT> correspondence_cell;
  size_t number_valid_correspondences = 0;

  if (pcl::isSamePointType<PointSource, PointTarget> ())
  {
    for (std::vector<int>::const_iterator input_index = indices_->begin (); input_index != indices_->end (); ++input_index)
    {
      if (target_correspondences_lookup_table_.getCorrespondence(input_->points[*input_index], max_distance_squared, correspondence_cell))
      {
        pcl::Correspondence& correspondence = correspondences[number_valid_correspondences++];
        correspondence.index_query = *input_index;
        correspondence.index_match = correspondence_cell.closest_point_index;
        correspondence.distance = correspondence_cell.distance_squared_to_closest_point;
      }
    }
  }
  else
  {
    PointTarget pt;
    for (std::vector<int>::const_iterator input_index = indices_->begin (); input_index != indices_->end (); ++input_index)
    {
      copyPoint (input_->points[*input_index], pt);
      if (target_correspondences_lookup_table_.getCorrespondence(input_->points[*input_index], max_distance_squared, correspondence_cell))
      {
        pcl::Correspondence& correspondence = correspondences[number_valid_correspondences++];
        correspondence.index_query = *input_index;
        correspondence.index_match = correspondence_cell.closest_point_index;
        correspondence.distance = correspondence_cell.distance_squared_to_closest_point;
      }
    }
  }
  correspondences.resize (number_valid_correspondences);
  deinitCompute ();

  PCL_DEBUG ("[pcl::registration::CorrespondenceEstimationLookupTable::determineCorrespondences] Computed %u correspondences (%u using the LookupTable and %u using the search tree)\n",
             number_valid_correspondences, target_correspondences_lookup_table_.getNumberOfQueriesOnLookupTable(), target_correspondences_lookup_table_.getNumberOfQueriesOnSearchTree());
  target_correspondences_lookup_table_.resetNumberOfQueriesOnLookupTable();
  target_correspondences_lookup_table_.resetNumberOfQueriesOnSearchTree();
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar, typename DistanceT> void
pcl::registration::CorrespondenceEstimationLookupTable<PointSource, PointTarget, Scalar, DistanceT>::determineReciprocalCorrespondences (
    pcl::Correspondences &correspondences, double max_distance)
{
  if (!initCompute ())
    return;

  if (!initComputeReciprocal())
    return;

  double max_distance_squared = max_distance * max_distance;
  correspondences.resize (indices_->size());
  pcl::registration::CorrespondenceLookupTableCell<DistanceT> correspondence_cell;
  pcl::registration::CorrespondenceLookupTableCell<DistanceT> reciprocal_correspondence_cell;
  unsigned int number_valid_correspondences = 0;

  if (pcl::isSamePointType<PointSource, PointTarget> ())
  {
    for (std::vector<int>::const_iterator input_index = indices_->begin (); input_index != indices_->end (); ++input_index)
    {
      if (target_correspondences_lookup_table_.getCorrespondence(input_->points[*input_index], max_distance_squared, correspondence_cell) &&
          source_correspondences_lookup_table_.getCorrespondence(target_->points[correspondence_cell.closest_point_index], max_distance_squared, reciprocal_correspondence_cell) &&
          *input_index == reciprocal_correspondence_cell.closest_point_index)
      {
        pcl::Correspondence& correspondence = correspondences[number_valid_correspondences++];
        correspondence.index_query = *input_index;
        correspondence.index_match = correspondence_cell.closest_point_index;
        correspondence.distance = correspondence_cell.distance_squared_to_closest_point;
      }
    }
  }
  else
  {
    PointTarget pt_src;
    for (std::vector<int>::const_iterator input_index = indices_->begin (); input_index != indices_->end (); ++input_index)
    {
      pcl::copyPoint (input_->points[*input_index], pt_src);
      if (target_correspondences_lookup_table_.getCorrespondence(input_->points[*input_index], max_distance_squared, correspondence_cell))
      {
        PointSource pt_tgt;
        pcl::copyPoint (target_->points[correspondence_cell.closest_point_index], pt_tgt);
        if (source_correspondences_lookup_table_.getCorrespondence(target_->points[correspondence_cell.closest_point_index], max_distance_squared, reciprocal_correspondence_cell) &&
            *input_index == reciprocal_correspondence_cell.closest_point_index)
        {
         pcl::Correspondence& correspondence = correspondences[number_valid_correspondences++];
          correspondence.index_query = *input_index;
          correspondence.index_match = correspondence_cell.closest_point_index;
          correspondence.distance = correspondence_cell.distance_squared_to_closest_point;
        }
      }
    }
  }
  correspondences.resize (number_valid_correspondences);
  deinitCompute ();

  PCL_DEBUG ("[pcl::registration::CorrespondenceEstimationLookupTable::determineReciprocalCorrespondences] Computed %u correspondences (%u using the LookupTable, %u using the search tree, %u using the reciprocal LookupTable and %u using the reciprocal search tree)\n",
               number_valid_correspondences, target_correspondences_lookup_table_.getNumberOfQueriesOnLookupTable(), target_correspondences_lookup_table_.getNumberOfQueriesOnSearchTree(),
               source_correspondences_lookup_table_.getNumberOfQueriesOnLookupTable(), source_correspondences_lookup_table_.getNumberOfQueriesOnSearchTree());
  target_correspondences_lookup_table_.resetNumberOfQueriesOnLookupTable();
  target_correspondences_lookup_table_.resetNumberOfQueriesOnSearchTree();
  source_correspondences_lookup_table_.resetNumberOfQueriesOnLookupTable();
  source_correspondences_lookup_table_.resetNumberOfQueriesOnSearchTree();
}

#endif /* PCL_REGISTRATION_IMPL_CORRESPONDENCE_ESTIMATION_LOOKUP_TABLE_H_ */
