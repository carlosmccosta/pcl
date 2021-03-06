/*
 * Software License Agreement (BSD License)
 *
 * Point Cloud Library (PCL) - www.pointclouds.org
 * Copyright (c) 2012-, Open Perception, Inc.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *  * Neither the name of the copyright holder(s) nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#pragma once

#include <pcl/filters/filter.h>
#include <memory>

namespace pcl
{
  /** \brief @b Scale Scales a point cloud by a factor
    *
    * \author Carlos M. Costa
    * \ingroup filters
    */
  template <typename PointT>
  class Scale: public Filter<PointT>
  {
    typedef typename Filter<PointT>::PointCloud PointCloud;

    using Filter<PointT>::filter_name_;
    using Filter<PointT>::input_;
    using Filter<PointT>::indices_;
    using Filter<PointT>::getClassName;

    public:
      using Ptr = std::shared_ptr<Scale<PointT> >;
      using ConstPtr = std::shared_ptr<const Scale<PointT> >;

      /** \brief Empty constructor. */
      Scale () :
        scale_ (0.001),
        transform_normals_ (true)
      {
        filter_name_ = "Scale";
      }

      /** \brief Destructor. */
      virtual ~Scale ()
      {
      }

      /** \brief Set the scale factor.
        * \param scale Scale factor
        */
      inline void
      setScale (double scale)
      { 
        scale_ = scale;
      }

      /** \brief Get the scale factor.
        * \param scale Scale factor
        */
      inline double
      getScale ()
      {
        return scale_;
      }

      /** \brief Set the transform normals flag.
        * \param transform_normals Transform normals flag
        */
      inline void
      setTransformNormals (bool transform_normals)
      {
        transform_normals_ = transform_normals;
      }

      /** \brief Get the transform normals flag.
        * \param scale Transform normals flag
        */
      inline bool
      getTransformNormals ()
      {
        return transform_normals_;
      }

    protected:
      /** \brief Scale factor. */
      double scale_;

      /** \brief Transform normals flag. */
      bool transform_normals_;

      /** \brief Scales a point cloud
        * \param[out] output the resultant point cloud
        */
      void
      applyFilter (PointCloud &output);
  };
}

#ifdef PCL_NO_PRECOMPILE
#include <pcl/filters/impl/scale.hpp>
#endif
