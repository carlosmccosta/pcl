/*
* Software License Agreement (BSD License)
*
*  Point Cloud Library (PCL) - www.pointclouds.org
*  Copyright (c) 2010-2011, Willow Garage, Inc.
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

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <pcl/point_types.h>
#include <pcl/registration/correspondence_estimation_lookup_table.h>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceLookupTable, ComputeLookupTableBounds)
{
  pcl::registration::CorrespondenceLookupTable<pcl::PointXYZ, float> lut;
  pcl::PointCloud<pcl::PointXYZ> pointcloud;
  ASSERT_FALSE (lut.computeLookupTableBounds(pointcloud));
  pointcloud.push_back (pcl::PointXYZ ( 0, 0, 0));
  pointcloud.push_back (pcl::PointXYZ (-1, 2,-3));
  pointcloud.push_back (pcl::PointXYZ ( 4,-5, 6));
  pointcloud.push_back (pcl::PointXYZ (-7, 8, 9));
  ASSERT_TRUE (lut.computeLookupTableBounds (pointcloud));

  Eigen::Vector3f margin = lut.getLookupTableMargin ();
  Eigen::Vector3f min_bounds = lut.getMinimumBounds ();
  Eigen::Vector3f max_bounds = lut.getMaximumBounds ();
  float half_cell_resolution = lut.getCellResolution () * 0.5;

  EXPECT_FLOAT_EQ ((float)(-7.0 - margin (0) - half_cell_resolution), (float)min_bounds (0));
  EXPECT_FLOAT_EQ ((float)(-5.0 - margin (1) - half_cell_resolution), (float)min_bounds (1));
  EXPECT_FLOAT_EQ ((float)(-3.0 - margin (2) - half_cell_resolution), (float)min_bounds (2));
  EXPECT_FLOAT_EQ ((float)( 4.0 + margin (0) + half_cell_resolution), (float)max_bounds (0));
  EXPECT_FLOAT_EQ ((float)( 8.0 + margin (1) + half_cell_resolution), (float)max_bounds (1));
  EXPECT_FLOAT_EQ ((float)( 9.0 + margin (2) + half_cell_resolution), (float)max_bounds (2));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceLookupTable, ComputeCorrespondenceCellIndex)
{
  pcl::registration::CorrespondenceLookupTable<pcl::PointXYZ, float> lut;
  lut.setCellResolution (1.0);
  lut.setLookupTableMargin (Eigen::Vector3f (0,0,0));
  lut.setUseSearchTreeWhenQueryPointIsOutsideLookupTable (false);
  typename pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud (new pcl::PointCloud<pcl::PointXYZ> ());
  pointcloud->push_back (pcl::PointXYZ (0,0,0));
  pointcloud->push_back (pcl::PointXYZ (2,3,4));
  ASSERT_TRUE (lut.initLookupTable (pointcloud));

  size_t correspondence_index;
  Eigen::Vector3i correspondence_index_components;
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1,-1,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1, 0,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1, 4,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1,-1, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1, 0, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1, 4, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1,-1, 5), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1, 0, 5), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (-1, 4, 5), correspondence_index, correspondence_index_components));

  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0,-1,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0, 0,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0, 4,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0,-1, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0, 4, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0,-1, 5), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0, 0, 5), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 0, 4, 5), correspondence_index, correspondence_index_components));

  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3,-1,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3, 0,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3, 4,-1), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3,-1, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3, 0, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3, 4, 0), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3,-1, 5), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3, 0, 5), correspondence_index, correspondence_index_components));
  EXPECT_FALSE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ ( 3, 4, 5), correspondence_index, correspondence_index_components));

  EXPECT_TRUE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (0,0,0), correspondence_index, correspondence_index_components));
  EXPECT_EQ (0, correspondence_index);
  EXPECT_EQ (Eigen::Vector3i (0,0,0), correspondence_index_components);

  EXPECT_TRUE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (1,1,1), correspondence_index, correspondence_index_components));
  EXPECT_EQ (1+3+3*4, correspondence_index);
  EXPECT_EQ (Eigen::Vector3i (1,1,1), correspondence_index_components);

  EXPECT_TRUE (lut.computeCorrespondenceCellIndex (pcl::PointXYZ (2,3,4), correspondence_index, correspondence_index_components));
  EXPECT_EQ (2+3*3+4*3*4, correspondence_index);
  EXPECT_EQ (Eigen::Vector3i (2,3,4), correspondence_index_components);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceLookupTable, ComputeCorrespondenceCellCentroid)
{
  pcl::registration::CorrespondenceLookupTable<pcl::PointXYZ, float> lut;
  lut.setCellResolution (1.0);
  lut.setLookupTableMargin (Eigen::Vector3f (0,0,0));
  typename pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud (new pcl::PointCloud<pcl::PointXYZ> ());
  pointcloud->push_back (pcl::PointXYZ (0,0,0));
  pointcloud->push_back (pcl::PointXYZ (2,3,4));
  ASSERT_TRUE (lut.initLookupTable (pointcloud));

  Eigen::Vector3f centroid = lut.computeCorrespondenceCellCentroid (Eigen::Vector3i (0,0,0));
  EXPECT_FLOAT_EQ (0, (float)centroid (0));
  EXPECT_FLOAT_EQ (0, (float)centroid (1));
  EXPECT_FLOAT_EQ (0, (float)centroid (2));

  centroid = lut.computeCorrespondenceCellCentroid (Eigen::Vector3i (1,0,0));
  EXPECT_FLOAT_EQ (1, (float)centroid (0));
  EXPECT_FLOAT_EQ (0, (float)centroid (1));
  EXPECT_FLOAT_EQ (0, (float)centroid (2));

  centroid = lut.computeCorrespondenceCellCentroid (Eigen::Vector3i (0,1,0));
  EXPECT_FLOAT_EQ (0, (float)centroid (0));
  EXPECT_FLOAT_EQ (1, (float)centroid (1));
  EXPECT_FLOAT_EQ (0, (float)centroid (2));

  centroid = lut.computeCorrespondenceCellCentroid (Eigen::Vector3i (0,0,1));
  EXPECT_FLOAT_EQ (0, (float)centroid (0));
  EXPECT_FLOAT_EQ (0, (float)centroid (1));
  EXPECT_FLOAT_EQ (1, (float)centroid (2));

  centroid = lut.computeCorrespondenceCellCentroid (Eigen::Vector3i (1,2,3));
  EXPECT_FLOAT_EQ (1, (float)centroid (0));
  EXPECT_FLOAT_EQ (2, (float)centroid (1));
  EXPECT_FLOAT_EQ (3, (float)centroid (2));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceLookupTable, GetCorrespondence)
{
  pcl::registration::CorrespondenceLookupTable<pcl::PointXYZ, float> lut;
  lut.setCellResolution (1.0);
  lut.setComputeDistanceFromQueryPointToClosestPoint (false);
  lut.setLookupTableMargin (Eigen::Vector3f (0,0,0));
  lut.setUseSearchTreeWhenQueryPointIsOutsideLookupTable (false);
  typename pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud (new pcl::PointCloud<pcl::PointXYZ> ());
  pointcloud->push_back (pcl::PointXYZ (0,0,0));
  pointcloud->push_back (pcl::PointXYZ (2,3,4));
  ASSERT_TRUE (lut.initLookupTable (pointcloud));

  pcl::registration::CorrespondenceLookupTableCell<float> correspondance;
  EXPECT_FALSE (lut.getCorrespondence (pcl::PointXYZ (-1, 0, 0), 4, correspondance));
  EXPECT_FALSE (lut.getCorrespondence (pcl::PointXYZ ( 0,-1, 0), 4, correspondance));
  EXPECT_FALSE (lut.getCorrespondence (pcl::PointXYZ ( 0, 0,-1), 4, correspondance));

  lut.setUseSearchTreeWhenQueryPointIsOutsideLookupTable (true);
  EXPECT_FALSE (lut.getCorrespondence (pcl::PointXYZ (-1, 0, 0), 0.99, correspondance));

  EXPECT_TRUE (lut.getCorrespondence (pcl::PointXYZ (-1, -2, -3), 14, correspondance));
  EXPECT_FLOAT_EQ (14, correspondance.distance_squared_to_closest_point);
  EXPECT_FLOAT_EQ (0, correspondance.closest_point_index);

  EXPECT_TRUE (lut.getCorrespondence (pcl::PointXYZ (2.1, 3.1, 4.1), 2, correspondance));
  EXPECT_FALSE(correspondance.distance_squared_to_closest_point != 0);
  EXPECT_FLOAT_EQ (1, correspondance.closest_point_index);

  lut.setComputeDistanceFromQueryPointToClosestPoint (true);
  EXPECT_TRUE (lut.getCorrespondence (pcl::PointXYZ (2, 3, 4), 2, correspondance));
  EXPECT_FLOAT_EQ (0, correspondance.distance_squared_to_closest_point);
  EXPECT_FLOAT_EQ (1, correspondance.closest_point_index);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceLookupTable, AddPointsToLookupTable)
{
  ASSERT_TRUE (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceLookupTable, InitLookupTable)
{
  ASSERT_TRUE (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceLookupTable, InitLookupTableUsingEuclideanDistanceTransform)
{
  ASSERT_TRUE (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceEstimationLookupTable, DetermineCorrespondencesKdtree)
{
  ASSERT_TRUE (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST (CorrespondenceEstimationLookupTable, DetermineCorrespondencesDistanceTransform)
{
  ASSERT_TRUE (true);
}

/* ---[ */
int
main (int argc, char** argv)
{
  testing::InitGoogleTest (&argc, argv);
  return (RUN_ALL_TESTS ());
}
/* ]--- */
