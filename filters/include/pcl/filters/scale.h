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

#ifndef PCL_FILTERS_SCALE_H_
#define PCL_FILTERS_SCALE_H_

#include <pcl/filters/filter.h>

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
      typedef boost::shared_ptr<Scale<PointT> > Ptr;
      typedef boost::shared_ptr<const Scale<PointT> > ConstPtr;

      /** \brief Empty constructor. */
      Scale () :
        scale_ (0.001)
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

      /** \brief Set the scale factor.
        * \param scale Scale factor
        */
      inline double
      getScale ()
      {
        return scale_;
      }

    protected:
      /** \brief Scale factor. */
      double scale_;

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

#endif  //#ifndef PCL_FILTERS_SCALE_H_

