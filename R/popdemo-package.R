################################################################################
#' Provides tools for demographic modelling using projection matrices
#'
#' @description
#' \code{popdemo} provides tools for modelling populations and demography using 
#' matrix projection models (MPMs), with deterministic and stochastic model 
#' implementations. These tools include population projection, indices of short- 
#' and long-term population size and growth, perturbation analysis, convergence 
#' to stability or stationarity, and diagnostic and manipulation tools. This includes:
#' \describe{
#'  \item{POPULATION PROJECTION}{
#'        \code{popdemo} provides a simple means of projecting and plotting PPM models.
#'        \code{\link{project}} provides a means to project and plot population 
#'        dynamics of both deterministic and stochastic models. Many methods are available
#'        for working with population projections: see \code{\link{Projection-class}} and 
#'        \code{\link{Projection-plots}}
#'  }
#'  \item{ASYMPTOTIC DYNAMICS}{
#'        The \code{\link{eigs}} function provides a simple means to calculate asymptotic 
#'        population dynamics using matrix eigenvalues.
#'  }
#'  \item{TRANSIENT DYNAMICS}{
#'        There are functions for calculating transient dynamics at various points of 
#'        the population projection. \code{\link{reac}} measures immediate transient 
#'        density of a population (within the first time step). \code{\link{maxamp}}, 
#'        \code{\link{maxatt}} are near-term indices that measure the largest and 
#'        smallest transient dynamics a population may exhibit overall, respectively. 
#'        \code{\link{inertia}} measures asymptotic population density relative to stable
#'        state, and has many perturbation methods in the package (see below). All transient 
#'        indices can be calculated using specific population structures, as well as bounds on
#'        population size.
#'  }
#'  \item{PERTURBATION ANALYSIS}{
#'        Methods for linear perturbation (sensitivity and elasticity) analysis of asymptotic dynamics 
#'        are available through the \code{\link{sens}}, \code{\link{tfs_lambda}} and \code{\link{tfsm_lambda}} 
#'        functions. Elasticity analysis is also available using the \code{\link{elas}} function. 
#'        Sensitivity analysis of transient dynamics is available using the \code{\link{tfs_inertia}} 
#'        and \code{\link{tfsm_inertia}} functions.
#'        Methods for nonlinear perturbation (transfer function) analysis of asymptotic 
#'        dynamics is achieved using \code{\link{tfa_lambda}} and \code{\link{tfam_lambda}}, 
#'        whilst transfer function analysis of transient dynamics is available with 
#'        \code{\link{tfa_inertia}} and \code{\link{tfam_inertia}}. These all have 
#'        associated plotting methods linked to them: see \code{\link{plot.tfa}} and 
#'        \code{\link{plot.tfam}}).
#'  }
#'  \item{MODEL CONVERGENCE}{
#'        Information on the convergence of populations to stable state can be useful, and
#'        \code{popdemo} provides several means of analysing convergence.
#'        \code{\link{dr}} measures the damping ratio, and there are several distance measures available
#'        (see \code{\link{KeyfitzD}}, \code{\link{projectionD}} and 
#'        \code{\link{CohenD}}).  There is also a means of calculating convergence time
#'        through simulation: \code{\link{convt}}.
#'  }
#'  \item{DIAGNOSTIC TOOLS}{
#'        \code{\link{isPrimitive}}, \code{\link{isIrreducible}} and
#'        \code{\link{isErgodic}} facilitate diagnosis of matrix properties
#'        pertaining to ergodicity.
#'  }
#'  \item{UTILITIES}{
#'        \code{\link{Matlab2R}} allows coding of matrices in a Matlab style, which
#'        also facilitates import of multiple matrices simultaneously if comma-seperated
#'        files are used to import dataframes. Its analogue, \code{\link{R2Matlab}},
#'        converts \R matrices to Matlab-style strings, for easier export.
#'  }
#' }
#'
"_PACKAGE"
