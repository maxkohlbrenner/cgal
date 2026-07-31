// Copyright (c) 2025  TU Berlin
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Max Kohlbrenner

#ifndef CGAL_MAXIMAL_EMPTY_SPHERES_MAXIMAL_EMPTY_SPHERES_H
#define CGAL_MAXIMAL_EMPTY_SPHERES_MAXIMAL_EMPTY_SPHERES_H

#include <CGAL/license/Maximal_empty_spheres.h>

#include <CGAL/Maximal_empty_spheres/internal/Lie_geometry.h>
#include <CGAL/Maximal_empty_spheres/internal/rotation_from_to.h>

#include <CGAL/Timer.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Dimension.h>

#include <Eigen/Core>

#include <fstream>
#include <iostream>
#include <vector>

namespace CGAL {

#ifndef DOXYGEN_RUNNING

template<typename Dimension>
void maximal_empty_spheres(const Eigen::MatrixXd& G,
                           Eigen::MatrixXd& result,
                           Eigen::MatrixXi* contact_indices = nullptr,
                           double atol = 1e-8,
                           int debug_level = 4)
{
  constexpr int D = Dimension::value;
  using namespace CGAL::Maximal_empty_spheres::internal;

  using Kernel = CGAL::Epick_d< CGAL::Dimension_tag<D+2> >;
  using Vertex = CGAL::Triangulation_vertex<Kernel, std::ptrdiff_t>;
  using DS_full_cell = CGAL::Triangulation_ds_full_cell<void,CGAL::TDS_full_cell_mirror_storage_policy>;
  using Full_cell = CGAL::Triangulation_full_cell<Kernel, int, DS_full_cell>;
  using TDS = CGAL::Triangulation_data_structure<CGAL::Dimension_tag<D+2>, Vertex, Full_cell>;
  using Triangulation = CGAL::Triangulation<Kernel,TDS>;

  CGAL_precondition(D == G.cols()-1 && "Dimension does not match the number of columns in G");

  bool full_simplices_only = true;
  int max_plains = 10;

  if (debug_level > 0) {
    std::cout << "Loaded SDF values of shape (" << G.rows() << ", " << G.cols() << ")"  << std::endl;
  }

  Eigen::MatrixXd SE = G;
  SE.col(D) = G.col(D).array().abs();

  // std::cout << "SE: " << std::endl;
  // std::cout << SE << std::endl;

  if (debug_level > 3) {
    std::cout << "Convert to Lie representation" << std::endl;
  }

  Eigen::MatrixXd H;
  lie_ip_matrix(D,H);

  Eigen::MatrixXd SL;
  spheres_to_lie(SE,SL);
  if (debug_level > 1) {
    std::cout << "Lie Repr: " << std::endl << SL << std::endl;
  }

  bool add_zerorad=false;
  if (add_zerorad){
    // add an additional sphere of zero radius that is completely hidden in another, non-zero sphere. 
      std::cout << "Addint a hidden sphere" << std::endl;
      Eigen::MatrixXd SLe(SL.rows()+1,SL.cols());
      SLe.block(0,0,SL.rows(),SL.cols()) = SL;
      // Eigen::RowVectorXd n_hsn = Eigen::RowVectorXd::Zero(NC.cols());
      SLe.row(SL.rows())       = SL.row(0);
      SLe(SL.rows(),SL.cols()-1)=0.;
      SL=SLe;
  }


  Eigen::MatrixXd NC = SL * H;
  // std::cout << "... on quadric? ";
  // std::cout << std::boolalpha << ( ((NC.array() * SL.array()).rowwise().sum()).matrix().norm() < atol) << std::endl;

  std::cout << "NC: " << NC << std::endl;

  bool rc_hs=false;
  if (rc_hs){
      std::cout << "Integrating the last component is negative halfspace" << std::endl;
      Eigen::MatrixXd NCe(NC.rows()+1,NC.cols());
      NCe.block(0,0,NC.rows(),NC.cols()) = NC;
      Eigen::RowVectorXd n_hsn = Eigen::RowVectorXd::Zero(NC.cols());
      n_hsn(SL.cols()-2)=-1;
      NCe.row(SL.rows())       = n_hsn;
      NC = NCe;
  }


  // std::cout << "NC: " << std::endl << NC << std::endl;

  Eigen::VectorXd s0;
  if (true){
      s0  = Eigen::VectorXd::Zero(D+3);
      s0(D+2) = -1.;
      s0(D)   = -1.;
  } else {
    s0 = NC.colwise().mean().transpose();
  }
  Eigen::VectorXd dst = Eigen::VectorXd::Zero(D+3);
  dst[D+2] = 1.;
  Eigen::MatrixXd R = rotation_from_to(s0, dst);

  // std::cout << "Rot: " << R << std::endl;

  if (debug_level > 0) {
    std::cout << std::boolalpha << "All normals in the same halfspace? : " <<  (((NC*s0).transpose()).array() >= atol).all() << std::endl;
  }

  Eigen::MatrixXd NCR = NC * R.transpose(); // rotate
  if (debug_level > 0) {
    std::cout << std::boolalpha << "Last component greater than zero?  : "  << (NCR.col(D+2).array()    >= atol).all() << std::endl;
  }
  // std::cout << "NCR:" << NCR << std::endl;

  NCR.array().colwise() /= NCR.array().col(D+2);
  Eigen::MatrixXd Np = NCR.block(0,0,NCR.rows(),D+2);

  // ------------------------ CONVEX HULL COMPUTATION -------------------------

  // --> convert to point list --
  std::vector<typename Triangulation::Point> points;
  points.reserve(int(Np.rows()));
  std::vector<typename Kernel::FT> p_;
  p_.reserve(D+2);
  for (long i=0; i<Np.rows(); ++i) {
    p_.clear();
    for (int d=0; d<D+2; ++d) {
      p_.push_back(Np(i,d));
    }
    points.push_back(typename Triangulation::Point(p_.begin(),p_.end()));
  }

  Timer timer;
  timer.start();

  Triangulation t(D+2);
  t.insert_and_index(points.begin(), points.end());

  if (debug_level > 0) {
    std::cout << "Convex hull of " << points.size() << " points computed in " << timer.time() << " seconds." << std::endl;
  }

  std::cout << "Current dimension: " << t.current_dimension() << std::endl;
  int current_dim = t.current_dimension();

  if (current_dim < D){
    std::cout << "[ERROR:] (current_dim < D) :TODO" << std::endl;
  } else if (current_dim == D){
    // problem of Appolonius, more configs possible? 
    std::cout << "[Not Implemented:] (current_dim == D) :TODO" << std::endl;
  } else if (current_dim == D+1){
    // WDT

    // get boundary facets of the CH and give them a unique index
    std::vector<typename Triangulation::Full_cell_handle> infinite_cells;
    for(auto it = t.full_cells_begin(); it != t.full_cells_end(); ++it) {
      if (t.is_infinite(it)){
          infinite_cells.push_back(it);
      }
    }
    int ci = 0;
    for(auto ch : infinite_cells) {
      ch->data() = ci++;
    }

    // check planes
    int n_check_planes = (NC.rows()<=max_plains)? NC.rows(): max_plains;
    Eigen::MatrixXd NC_(n_check_planes,D+3);
    for (int i=0; i<n_check_planes; i++){
        NC_.row(i) = NC.row((i*20747) % NC.rows());
    } 

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    Eigen::MatrixXd simplex(current_dim,D+3);
    Eigen::RowVectorXi simplex_inds(current_dim);
    Eigen::MatrixXd K_l(2,D+3);
    Eigen::Vector<bool,Eigen::Dynamic>  inv;

    std::vector<Eigen::RowVectorXd> solutions_;
    std::vector<Eigen::RowVectorXi> contact_indices_;

    for(auto ch : infinite_cells) {

        // extract simplex
        int si=0;
        for (int i=0; i<current_dim+1; ++i) {
            std::cout << "  " << i << " " << std::endl;
            if(ch->vertex(i) != t.infinite_vertex()) {
                typename Triangulation::Vertex_handle vh = ch->vertex(i);
                int vid = vh->data();
                Eigen::RowVectorXd nci = NC.row(vid);
                simplex.row(si)  = nci;
                simplex_inds(si) = vid;
                si++;
            } else std::cout << " --INF--" << std::endl;
        }
        std::cout << "Simplex: " << std::endl << simplex << std::endl;
        std::cout << "(inds) : " << std::endl << simplex_inds << std::endl;
    
        // span the normal space and flip into the cone
        svd.compute(simplex, Eigen::ComputeFullV);

        std::cout << "svalues: " << svd.singularValues().transpose() << std::endl;

        K_l = svd.matrixV().block(0,D+1,D+3,2).transpose();
        inv = ((K_l * NC_.transpose()).rowwise().maxCoeff().array() >= atol);
        for (int i=0; i<K_l.rows(); ++i) {
          if (inv(i)) {
            K_l.row(i) *= -1;
          }
        }
        std::cout << "In Cone? " << std::endl;
        std::cout << K_l * NC_.transpose() << std::endl;
        std::cout << std::endl;

        double ls[2]; // 1,l2;
        line_quadric_intersection(K_l.row(0).transpose(), K_l.row(1).transpose(), H, ls[0],ls[1]);
        std::cout << "ls: " << ls[0] << ", " << ls[1] << std::endl;

        for (int li=0; li<2; ++li) {
            // if ((0.-atol<=ls[li]) && (ls[li]<=1.+atol)) {
            if (true) {
                std::cout << "cand: " << ls[li] << std::endl;
                Eigen::RowVectorXd s_ = (1-ls[li])*K_l.row(0)+ls[li]*K_l.row(1);
                std::cout << s_ << std::endl;
                if ((s_(D+2) < 0.) && (s_(D+1) >= 0) && (fabs(s_(D+2)) >= atol)) {
                // if (true) {
                    std::cout << "--> Sol" << std::endl;
                    std::cout << s_ << std::endl;
                    solutions_.push_back(s_);
                    if (contact_indices) {
                      contact_indices_.emplace_back(simplex_inds);
                    }
                }
            }
        }
    }

    Eigen::MatrixXd solutions(solutions_.size(),D+3);
    for (std::size_t i=0; i<solutions_.size(); ++i) {
      solutions.row(i) = solutions_[i];
    }

    std::cout << "Solutions in Cone? " << std::endl;
    std::cout << solutions * NC_.transpose() << std::endl;


    lie_to_spheres(solutions, result);

    if(debug_level > 0) {
      std::cout << "Solutions.size: " << solutions.rows() << ", " << solutions.cols() << std::endl;
    }

    if (contact_indices) {
      contact_indices->resize(solutions_.size(),D+1);
      for (std::size_t i=0; i<solutions_.size(); ++i) contact_indices->block(i,0,1,D+1) = contact_indices_[i];
    }

    std::cout << "[Not Implemented:] (current_dim == D+1) :TODO" << std::endl;

  } else {
    // standard case

    // --> recover the infinite cells (contain the convex hull)
    std::vector<typename Triangulation::Full_cell_handle> infinite_cells;
    for(auto it = t.full_cells_begin(); it != t.full_cells_end(); ++it) {
      // bool infcell = t.is_infinite(it);
      // if( (t.is_infinite(it) && (current_dim == D+2) ) ||
      //    !t.is_infinite(it) && (current_dim == D+1) ) {
      if (t.is_infinite(it)){
          infinite_cells.push_back(it);
      }
    }

    // --> get a unique index per cell
    int ci = 0;
    for(auto ch : infinite_cells) {
    ch->data() = ci++;
    }

    if (debug_level > 0) {
    std::cout << "Boundary of CH has " << infinite_cells.size() << " facets" << std::endl;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    Eigen::MatrixXd Ks(infinite_cells.size(),D+3);
    Eigen::Vector<bool,Eigen::Dynamic> full_simplices(infinite_cells.size());

    Eigen::MatrixXd simplex(current_dim,D+3);
    Eigen::RowVectorXi simplex_inds(current_dim);
    int ki = 0;
    int nafs = 0;

    for(auto ch : infinite_cells) {
    std::cout << "ch: " << *ch << std::endl;
    int si=0;
    for (int i=0; i<current_dim+1; ++i) {
        std::cout << "  " << i << " " << std::endl;
        if(ch->vertex(i) != t.infinite_vertex()) {
            typename Triangulation::Vertex_handle vh = ch->vertex(i);
            // std::cout << "        vh: " << *vh << std::endl;
            int vid = vh->data();
            // std::cout << "        vid: " << vid << std::endl;
            Eigen::RowVectorXd nci = NC.row(vid);
            simplex.row(si)  = nci;
            simplex_inds(si) = vid;
            si++;
        } else std::cout << " --INF--" << std::endl;
    }

    std::cout << "Simplex inds: " << simplex_inds << std::endl;

    if (si == D+1){
        // WDT case
        std::cout << "Simplex: " << std::endl << simplex << std::endl;
        std::cout << "(inds) : " << std::endl << simplex_inds << std::endl;
        svd.compute(simplex, Eigen::ComputeFullV);
        std::cout << svd.matrixV().cols() << std::endl;

        double ls[2]; // 1,l2;
        line_quadric_intersection(svd.matrixV().col(D+1), svd.matrixV().col(D+2), H, ls[0],ls[1]);
        std::cout << "ls: " << ls[0] << ", " << ls[1] << std::endl;
        for (int li=0; li<2; ++li) {
            Eigen::RowVectorXd s_ = (1-ls[li])*svd.matrixV().col(D+1)+ls[li]*svd.matrixV().col(D+2);
            if ((s_(D+2) < 0.) && (s_(D+1) >= 0) && (fabs(s_(D+2)) >= atol)) {
                std::cout << s_ << std::endl;
            }
        }


    } else if (si < D+2){
        // std::cout << "NFS BOUNDARY SIMPLEX" << std::endl;
        full_simplices(ki) = false;
        ++nafs;
        ++ki;
        std::cout << "(inds) : " << std::endl << simplex_inds << std::endl;
    } else {
        std::cout << "Simplex: " << std::endl << simplex << std::endl;
        std::cout << "(inds) : " << std::endl << simplex_inds << std::endl;
        svd.compute(simplex, Eigen::ComputeFullV);
        Ks.row(ki)         = svd.matrixV().col(D+2).transpose();
        full_simplices(ki) = svd.singularValues().array().abs().minCoeff() > atol;

        std::cout << "svalues: " << svd.singularValues().transpose() << std::endl;
        std::cout << "vertex:  " << Ks.row(ki) << std::endl;

        ++ki;

        if(debug_level > 0) {
          if (svd.singularValues().array().abs().minCoeff() <= atol) {
            ++nafs;
          }
        }
    }

    }

    if (debug_level > 0 && nafs > 0) {
      std::cout << "Warning: " << nafs << " simplices were not full simplices, i.e. they are not full-dimensional in the ambient space." << std::endl;
    }

    // std::cout << "Ks.shape: (" << Ks.rows() << ", " << Ks.cols() << ")" << std::endl;
    // std::cout << "Ks: " << std::endl << Ks << std::endl;
    // std::cout << "Ks * NC.T" << std::endl << Ks * NC.transpose() << std::endl;

    // @todo: reduced to 200 only for memory reasons, loop and increase again? do sth smarter?
    int n_check_planes = (NC.rows()<=max_plains)? NC.rows(): max_plains;
    if(debug_level > 0) {
    std::cout << "n_check_planes: " << n_check_planes << std::endl;
    std::cout << "NC.shape(): " << NC.rows() << ", " << NC.cols() << std::endl;
    }

    Eigen::Vector<bool,Eigen::Dynamic>  inv;
    if (true){
    Eigen::MatrixXd NC_(n_check_planes,D+3);
    for (int i=0; i<n_check_planes; i++){
        NC_.row(i) = NC.row((i*20747) % NC.rows());
    } 
      inv = ((Ks * NC_.transpose()).rowwise().maxCoeff().array() >= atol);
    } else {
      inv = ((Ks * NC.block(0,0,n_check_planes,D+3).transpose()).rowwise().maxCoeff().array() >= atol);
    }
    for (int i=0; i<Ks.rows(); ++i) {
      if (inv(i)) {
        Ks.row(i) *= -1;
      }
    }

    if(debug_level > 0) {
    std::cout << "... plane flips done" << std::endl;
    }

    if (debug_level > 3) {
    // costly and not possible for large sets of spheres
    std::cout << "Debugging, remove later (matrices too large), all Ks in the cone?: ";
    std::cout << ((NC * Ks.transpose()).rowwise().maxCoeff().array() <= atol).array().all() << std::endl;

    std::cout << "Ks: " << std::endl;
    std::cout << Ks << std::endl;
    std::cout << std::endl;
    std::cout << "valids: " << full_simplices.transpose() << std::endl;
    std::cout << std::endl;
    }

    bool cone_filter=false;
    int n_cone_filtered=0;
    std::vector<Eigen::RowVectorXd> solutions_;
    std::vector<Eigen::RowVectorXi> contact_indices_;
    for(auto ch : infinite_cells) {
    int ci = ch->data();

    for (int i=0; i<current_dim+1; ++i) {
      if(ch->vertex(i) != t.infinite_vertex()) {
          /*
        if(((ch->vertex(i)         != t.infinite_vertex()) && (current_dim == D+2) ) ||
             !t.is_infinite(ch->neighbor(i))               && (current_dim  < D+2) ){
            */
        int cj =  ch->neighbor(i)->data();
        if ((ci < cj) &&    // only treat each edge once
            (!full_simplices_only || (full_simplices(ci) && full_simplices(cj)))) {  // only consider edges between full simplices

          std::cout << "--- Edge (" << ci << " - " << cj << ") ---" << std::endl;
          std::cout << Ks.row(ci) << std::endl;
          std::cout << Ks.row(cj) << std::endl;

          double ls[2]; // 1,l2;
          line_quadric_intersection(Ks.row(ci), Ks.row(cj), H, ls[0],ls[1]);
          std::cout << "ls: " << ls[0] << ", " << ls[1] << std::endl;
          for (int li=0; li<2; ++li) {
            if ((0.-atol<=ls[li]) && (ls[li]<=1.+atol)) {
              std::cout << "POTENTIALY ADDING" << std::endl;
              Eigen::RowVectorXd s_ = (1-ls[li])*Ks.row(ci)+ls[li]*Ks.row(cj);

              if ((s_(D+2) < 0.) && (s_(D+1) >= 0) && (fabs(s_(D+2)) >= atol)) {
                bool add = true;
                if (cone_filter) {
                  if (((s_*NC.transpose()).array() > atol).any()) {
                    add = false;
                    ++n_cone_filtered;
                  }
                }

                std::cout << "add: " << add << std::endl;

                if (add) {
                  solutions_.push_back(s_);
                  if (contact_indices) {
                    Eigen::RowVectorXi st(D+1);
                    int ni=0;
                    for (int j=0; j<D+2; ++j) {
                      if (ch->vertex((i+1+j)%(D+3)) != t.infinite_vertex()) {
                        st(ni++) = ch->vertex((i+1+j)%(D+3))->data();
                      }
                    }
                    contact_indices_.emplace_back(st);
                  }
                }
              }
            }
          }
        }
      }
    }
    }

    if(debug_level > 0 && n_cone_filtered > 0) {
    std::cout << "(Cone filter!!!) " << std::endl;
    std::cout << " Filtered " << n_cone_filtered << " spheres that were not in the cone" << std::endl;
    }

    Eigen::MatrixXd solutions(solutions_.size(),D+3);
    for (std::size_t i=0; i<solutions_.size(); ++i) {
    solutions.row(i) = solutions_[i];
    }

    lie_to_spheres(solutions, result);

    if(debug_level > 0) {
    std::cout << "Solutions.size: " << solutions.rows() << ", " << solutions.cols() << std::endl;
    }

    if (contact_indices) {
    contact_indices->resize(solutions_.size(),D+1);
    for (std::size_t i=0; i<solutions_.size(); ++i) contact_indices->block(i,0,1,D+1) = contact_indices_[i];
    }

    }

    }

    template<typename SphereRange, typename OutputIterator, class DimensionTag>
    void maximal_empty_spheres(const SphereRange& /* input*/ ,
                           OutputIterator /* result*/,
                           DimensionTag /* tag */ )
    {
    CGAL_assertion_msg(false, "maximal_empty_spheres() called with unsupported dimension tag");
    }

    template<typename SphereRange, typename OutputIterator>
    void maximal_empty_spheres(const SphereRange& input,
                           OutputIterator result,
                           CGAL::Dimension_tag<3> /* tag */)
    {
    using Sphere_3 = typename SphereRange::value_type;
    using Kernel = typename Kernel_traits<Sphere_3>::Kernel;
    using Point_3 = typename Kernel::Point_3;

    Eigen::MatrixXd G(input.size(),4), res;
    for(auto it = input.begin(); it != input.end(); ++it) {
      G.row(it - input.begin()) = Eigen::RowVector4d(it->center().x(), it->center().y(), it->center().z(), sqrt(it->squared_radius()));
    }

    maximal_empty_spheres<Dimension_tag<3>>(G, res);
    for(int i=0; i<res.rows(); ++i) {
    Point_3 p(res(i,0), res(i,1), res(i,2));
    *result++ = Sphere_3(p, res(i,3)*res(i,3)); // res(i,3) is the radius
    }

#ifdef CGAL_MES_WRITE_SPHERES
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    // write soution in lie representation:  AF: I think this is not Lie representation, but rather the sphere representation
    std::ofstream pointfile("solution_lie_3.csv");
    pointfile << res.format(CSVFormat);
#endif
    }

    template<typename SphereRange, typename OutputIterator>
    void maximal_empty_spheres(const SphereRange& input,
                           OutputIterator result,
                           CGAL::Dimension_tag<2> /* tag */)
    {
    using Circle_2 = typename SphereRange::value_type;
    using Kernel = typename Kernel_traits<Circle_2>::Kernel;
    using Point_2 = typename Kernel::Point_2;

    Eigen::MatrixXd G(input.size(),3), res;
    for(auto it = input.begin(); it != input.end(); ++it) {
    G.row(it - input.begin()) = Eigen::RowVector3d(it->center().x(), it->center().y(), sqrt(it->squared_radius()));
    }

    maximal_empty_spheres<Dimension_tag<2>>(G, res);
    for(int i=0; i<res.rows(); ++i) {
    Point_2 p(res(i,0), res(i,1));
    *result++ = Circle_2(p, res(i,2)*res(i,2)); // res(i,2) is the radius
    }

#ifdef CGAL_MES_WRITE_SPHERES
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    // write soution in lie representation:
    std::ofstream pointfile("solution_lie_2.csv");
    pointfile << res.format(CSVFormat);
#endif
    }

#endif // DOXYGEN_RUNNING

    /*!
    * \ingroup PkgMaximalEmptySpheresFunctions
    *
    * \brief compute maximal empty spheres from a range of spheres.
    *
    * This function computes maximal empty spheres from the input range of spheres and stores the results in the output iterator.
    *
    * \tparam SphereRange A range of spheres to be processed. The value type may be a 2D circle or a 3D or dD sphere.
    * \tparam OutputIterator An output iterator to store the results which are the same type as the input spheres.
    *
    * \param input The input range of spheres.
    * \param result The output iterator to store the maximal empty spheres.
    */
    template <typename SphereRange, typename OutputIterator>
    void maximal_empty_spheres(const SphereRange& input, OutputIterator result)
    {
    using Sphere = typename SphereRange::value_type;

    constexpr int D = Ambient_dimension<Sphere>::value;
    maximal_empty_spheres(input, result, CGAL::Dimension_tag<D>());
    }

} // namespace CGAL

#endif // CGAL_MAXIMAL_EMPTY_SPHERES_MAXIMAL_EMPTY_SPHERES_H
