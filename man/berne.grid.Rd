\name{berne.grid}
\alias{berne.grid}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Berne -- very small extract of prediction grid
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The Berne Grid data set contains values of spatial covariates on the nodes of a 20 m grid. The data set is intended for spatial continouous predictions of soil properties modelled from the sampling locations in \code{\link{berne}} data set.
}
\usage{data("berne")}
\format{
  A data frame with 4594 observations on the following 228 variables.
  %% data.frame( paste0( "\\item{\\code{", names(dat), "}}{}" ) )
  \describe{
      \item{\code{id}}{node identifier number.}
      \item{\code{x}}{easting, Swiss grid in m, EPSG: 21781 (CH1903/LV03)}
      \item{\code{y}}{northing, Swiss grid in m, EPSG: 21781 (CH1903/LV03)}
 \item{\code{cl_mt_etap_pe}}{columns 4 to 228 contain environmental covariates representing soil forming factors. For more information see \code{Details} in \code{\link{berne}}.}
 \item{\code{cl_mt_etap_ro}}{}
 \item{\code{cl_mt_gh_1}}{}
 \item{\code{cl_mt_gh_10}}{}
 \item{\code{cl_mt_gh_11}}{}
 \item{\code{cl_mt_gh_12}}{}
 \item{\code{cl_mt_gh_2}}{}
 \item{\code{cl_mt_gh_3}}{}
 \item{\code{cl_mt_gh_4}}{}
 \item{\code{cl_mt_gh_5}}{}
 \item{\code{cl_mt_gh_6}}{}
 \item{\code{cl_mt_gh_7}}{}
 \item{\code{cl_mt_gh_8}}{}
 \item{\code{cl_mt_gh_9}}{}
 \item{\code{cl_mt_gh_y}}{}
 \item{\code{cl_mt_pet_pe}}{}
 \item{\code{cl_mt_pet_ro}}{}
 \item{\code{cl_mt_rr_1}}{}
 \item{\code{cl_mt_rr_10}}{}
 \item{\code{cl_mt_rr_11}}{}
 \item{\code{cl_mt_rr_12}}{}
 \item{\code{cl_mt_rr_2}}{}
 \item{\code{cl_mt_rr_3}}{}
 \item{\code{cl_mt_rr_4}}{}
 \item{\code{cl_mt_rr_5}}{}
 \item{\code{cl_mt_rr_6}}{}
 \item{\code{cl_mt_rr_7}}{}
 \item{\code{cl_mt_rr_8}}{}
 \item{\code{cl_mt_rr_9}}{}
 \item{\code{cl_mt_rr_y}}{}
 \item{\code{cl_mt_swb_pe}}{}
 \item{\code{cl_mt_swb_ro}}{}
 \item{\code{cl_mt_td_1}}{}
 \item{\code{cl_mt_td_10}}{}
 \item{\code{cl_mt_td_11}}{}
 \item{\code{cl_mt_td_12}}{}
 \item{\code{cl_mt_td_2}}{}
 \item{\code{cl_mt_tt_1}}{}
 \item{\code{cl_mt_tt_11}}{}
 \item{\code{cl_mt_tt_12}}{}
 \item{\code{cl_mt_tt_3}}{}
 \item{\code{cl_mt_tt_4}}{}
 \item{\code{cl_mt_tt_5}}{}
 \item{\code{cl_mt_tt_6}}{}
 \item{\code{cl_mt_tt_7}}{}
 \item{\code{cl_mt_tt_8}}{}
 \item{\code{cl_mt_tt_9}}{}
 \item{\code{cl_mt_tt_y}}{}
 \item{\code{ge_caco3}}{}
 \item{\code{ge_geo500h1id}}{}
 \item{\code{ge_geo500h3id}}{}
 \item{\code{ge_gt_ch_fil}}{}
 \item{\code{ge_lgm}}{}
 \item{\code{ge_vszone}}{}
 \item{\code{sl_nutr_fil}}{}
 \item{\code{sl_physio_neu}}{}
 \item{\code{sl_retention_fil}}{}
 \item{\code{sl_skelett_r_fil}}{}
 \item{\code{sl_wet_fil}}{}
 \item{\code{tr_be_gwn25_hdist}}{}
 \item{\code{tr_be_gwn25_vdist}}{}
 \item{\code{tr_be_twi2m_7s_tcilow}}{}
 \item{\code{tr_be_twi2m_s60_tcilow}}{}
 \item{\code{tr_ch_3_80_10}}{}
 \item{\code{tr_ch_3_80_10s}}{}
 \item{\code{tr_ch_3_80_20s}}{}
 \item{\code{tr_cindx10_25}}{}
 \item{\code{tr_cindx50_25}}{}
 \item{\code{tr_curv_all}}{}
 \item{\code{tr_curv_plan}}{}
 \item{\code{tr_curv_prof}}{}
 \item{\code{tr_enessk}}{}
 \item{\code{tr_es25}}{}
 \item{\code{tr_flowlength_up}}{}
 \item{\code{tr_global_rad_ch}}{}
 \item{\code{tr_lsf}}{}
 \item{\code{tr_mrrtf25}}{}
 \item{\code{tr_mrvbf25}}{}
 \item{\code{tr_ndom_veg2m_fm}}{}
 \item{\code{tr_nego}}{}
 \item{\code{tr_nnessk}}{}
 \item{\code{tr_ns25}}{}
 \item{\code{tr_ns25_145mn}}{}
 \item{\code{tr_ns25_145sd}}{}
 \item{\code{tr_ns25_75mn}}{}
 \item{\code{tr_ns25_75sd}}{}
 \item{\code{tr_poso}}{}
 \item{\code{tr_protindx}}{}
 \item{\code{tr_se_alti10m_c}}{}
 \item{\code{tr_se_alti25m_c}}{}
 \item{\code{tr_se_alti2m_fmean_10c}}{}
 \item{\code{tr_se_alti2m_fmean_25c}}{}
 \item{\code{tr_se_alti2m_fmean_50c}}{}
 \item{\code{tr_se_alti2m_fmean_5c}}{}
 \item{\code{tr_se_alti2m_std_10c}}{}
 \item{\code{tr_se_alti2m_std_25c}}{}
 \item{\code{tr_se_alti2m_std_50c}}{}
 \item{\code{tr_se_alti2m_std_5c}}{}
 \item{\code{tr_se_alti50m_c}}{}
 \item{\code{tr_se_alti6m_c}}{}
 \item{\code{tr_se_conv2m}}{}
 \item{\code{tr_se_curv10m}}{}
 \item{\code{tr_se_curv25m}}{}
 \item{\code{tr_se_curv2m}}{}
 \item{\code{tr_se_curv2m_s15}}{}
 \item{\code{tr_se_curv2m_s30}}{}
 \item{\code{tr_se_curv2m_s60}}{}
 \item{\code{tr_se_curv2m_s7}}{}
 \item{\code{tr_se_curv2m_std_10c}}{}
 \item{\code{tr_se_curv2m_std_25c}}{}
 \item{\code{tr_se_curv2m_std_50c}}{}
 \item{\code{tr_se_curv2m_std_5c}}{}
 \item{\code{tr_se_curv50m}}{}
 \item{\code{tr_se_curv6m}}{}
 \item{\code{tr_se_curvplan10m}}{}
 \item{\code{tr_se_curvplan25m}}{}
 \item{\code{tr_se_curvplan2m}}{}
 \item{\code{tr_se_curvplan2m_grass_17c}}{}
 \item{\code{tr_se_curvplan2m_grass_45c}}{}
 \item{\code{tr_se_curvplan2m_grass_9c}}{}
 \item{\code{tr_se_curvplan2m_s15}}{}
 \item{\code{tr_se_curvplan2m_s30}}{}
 \item{\code{tr_se_curvplan2m_s60}}{}
 \item{\code{tr_se_curvplan2m_s7}}{}
 \item{\code{tr_se_curvplan2m_std_10c}}{}
 \item{\code{tr_se_curvplan2m_std_25c}}{}
 \item{\code{tr_se_curvplan2m_std_50c}}{}
 \item{\code{tr_se_curvplan2m_std_5c}}{}
 \item{\code{tr_se_curvplan50m}}{}
 \item{\code{tr_se_curvplan6m}}{}
 \item{\code{tr_se_curvprof10m}}{}
 \item{\code{tr_se_curvprof25m}}{}
 \item{\code{tr_se_curvprof2m}}{}
 \item{\code{tr_se_curvprof2m_grass_17c}}{}
 \item{\code{tr_se_curvprof2m_grass_45c}}{}
 \item{\code{tr_se_curvprof2m_grass_9c}}{}
 \item{\code{tr_se_curvprof2m_s15}}{}
 \item{\code{tr_se_curvprof2m_s30}}{}
 \item{\code{tr_se_curvprof2m_s60}}{}
 \item{\code{tr_se_curvprof2m_s7}}{}
 \item{\code{tr_se_curvprof2m_std_10c}}{}
 \item{\code{tr_se_curvprof2m_std_25c}}{}
 \item{\code{tr_se_curvprof2m_std_50c}}{}
 \item{\code{tr_se_curvprof2m_std_5c}}{}
 \item{\code{tr_se_curvprof50m}}{}
 \item{\code{tr_se_curvprof6m}}{}
 \item{\code{tr_se_diss2m_10c}}{}
 \item{\code{tr_se_diss2m_25c}}{}
 \item{\code{tr_se_diss2m_50c}}{}
 \item{\code{tr_se_diss2m_5c}}{}
 \item{\code{tr_se_e_aspect10m}}{}
 \item{\code{tr_se_e_aspect25m}}{}
 \item{\code{tr_se_e_aspect2m}}{}
 \item{\code{tr_se_e_aspect2m_10c}}{}
 \item{\code{tr_se_e_aspect2m_25c}}{}
 \item{\code{tr_se_e_aspect2m_50c}}{}
 \item{\code{tr_se_e_aspect2m_5c}}{}
 \item{\code{tr_se_e_aspect2m_grass_17c}}{}
 \item{\code{tr_se_e_aspect2m_grass_45c}}{}
 \item{\code{tr_se_e_aspect2m_grass_9c}}{}
 \item{\code{tr_se_e_aspect50m}}{}
 \item{\code{tr_se_e_aspect6m}}{}
 \item{\code{tr_se_mrrtf2m}}{}
 \item{\code{tr_se_mrvbf2m}}{}
 \item{\code{tr_se_n_aspect10m}}{}
 \item{\code{tr_se_n_aspect25m}}{}
 \item{\code{tr_se_n_aspect2m}}{}
 \item{\code{tr_se_n_aspect2m_10c}}{}
 \item{\code{tr_se_n_aspect2m_25c}}{}
 \item{\code{tr_se_n_aspect2m_50c}}{}
 \item{\code{tr_se_n_aspect2m_5c}}{}
 \item{\code{tr_se_n_aspect2m_grass_17c}}{}
 \item{\code{tr_se_n_aspect2m_grass_45c}}{}
 \item{\code{tr_se_n_aspect2m_grass_9c}}{}
 \item{\code{tr_se_n_aspect50m}}{}
 \item{\code{tr_se_n_aspect6m}}{}
 \item{\code{tr_se_no2m_r500}}{}
 \item{\code{tr_se_po2m_r500}}{}
 \item{\code{tr_se_rough2m_10c}}{}
 \item{\code{tr_se_rough2m_25c}}{}
 \item{\code{tr_se_rough2m_50c}}{}
 \item{\code{tr_se_rough2m_5c}}{}
 \item{\code{tr_se_rough2m_rect3c}}{}
 \item{\code{tr_se_sar2m}}{}
 \item{\code{tr_se_sca2m}}{}
 \item{\code{tr_se_slope10m}}{}
 \item{\code{tr_se_slope25m}}{}
 \item{\code{tr_se_slope2m}}{}
 \item{\code{tr_se_slope2m_grass_17c}}{}
 \item{\code{tr_se_slope2m_grass_45c}}{}
 \item{\code{tr_se_slope2m_grass_9c}}{}
 \item{\code{tr_se_slope2m_s15}}{}
 \item{\code{tr_se_slope2m_s30}}{}
 \item{\code{tr_se_slope2m_s60}}{}
 \item{\code{tr_se_slope2m_s7}}{}
 \item{\code{tr_se_slope2m_std_10c}}{}
 \item{\code{tr_se_slope2m_std_25c}}{}
 \item{\code{tr_se_slope2m_std_50c}}{}
 \item{\code{tr_se_slope2m_std_5c}}{}
 \item{\code{tr_se_slope50m}}{}
 \item{\code{tr_se_slope6m}}{}
 \item{\code{tr_se_toposcale2m_r3_r50_i10s}}{}
 \item{\code{tr_se_tpi_2m_10c}}{}
 \item{\code{tr_se_tpi_2m_25c}}{}
 \item{\code{tr_se_tpi_2m_50c}}{}
 \item{\code{tr_se_tpi_2m_5c}}{}
 \item{\code{tr_se_tri2m_altern_3c}}{}
 \item{\code{tr_se_tsc10_2m}}{}
 \item{\code{tr_se_twi2m}}{}
 \item{\code{tr_se_twi2m_s15}}{}
 \item{\code{tr_se_twi2m_s30}}{}
 \item{\code{tr_se_twi2m_s60}}{}
 \item{\code{tr_se_twi2m_s7}}{}
 \item{\code{tr_se_vrm2m}}{}
 \item{\code{tr_se_vrm2m_r10c}}{}
 \item{\code{tr_slope25_l2g}}{}
 \item{\code{tr_terrtextur}}{}
 \item{\code{tr_tpi2000c}}{}
 \item{\code{tr_tpi5000c}}{}
 \item{\code{tr_tpi500c}}{}
 \item{\code{tr_tsc25_18}}{}
 \item{\code{tr_tsc25_40}}{}
 \item{\code{tr_twi2}}{}
 \item{\code{tr_twi_normal}}{}
 \item{\code{tr_vdcn25}}{}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
Due to CRAN file size restrictions the grid for spatial only shows a very small excerpt of the original study area.

The environmental covariates for prediction of soil properties from dataset \code{\link{berne}} were extracted at the nodes of a 20 m grid. For higher resolution geodata sets no averaging over the area of the 20x20 pixel was done. \code{Berne.grid} therefore has the same spatial support for each covariate as \code{\link{berne}}.

For more information on the environmental covariates see \code{\link{berne}}.

}
\references{
%%  ~~ possibly secondary sources and usages ~~
Nussbaum et al. (in prep). Comparison of statistical approaches for digital soil mapping.

}
\examples{

\dontrun{
data(berne.grid)

# plot spatial object
library(raster)

coordinates(berne.grid) <- ~x+y
proj4string(berne.grid) <- CRS("+init=epsg:21781")
gridded(berne.grid) <- TRUE

plot( raster(berne.grid, layer = "tr_se_mrrtf2m"))
}
}
\keyword{datasets}
