\name{berne}
\alias{berne}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Berne -- soil mapping case study
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The Berne dataset contains soil responses and a large set of explanatory covariates. The study area is located to the Northwest of the city of Berne and covers agricultural area.
Soil responses included are soil pH (4 depth intervals calculated from soil horizon), drainage classes (3 ordered classes) and presence of waterlogging characteristics down to a specified depth (binary response).

Covariates cover environmental conditions by representing climate, topography, parent material and soil.

}
\usage{data("berne")}
\format{
  A data frame with 1052 observations on the following 238 variables.
  %% data(berne); cat( paste0( "\\item{\\code{", names(berne), "}}{}\n" ))
  \describe{
      \item{\code{site_id_unique}}{ID of original profile sampling}
      \item{\code{x}}{easting, Swiss grid in m, EPSG: 21781 (CH1903/LV03)}
      \item{\code{y}}{northing, Swiss grid in m, EPSG: 21781 (CH1903/LV03)}
      \item{\code{dataset}}{Factor splitting dataset for \code{calibration} and independent \code{validation}. \code{validation} was assigned at random by using weights to ensure even spatial coverage of the agricultural area.}
      \item{\code{dclass}}{Drainage class, ordered Factor.}
      \item{\code{waterlog.30}}{Presence of waterlogging characteristics down to 30 cm (1: presence, 0: absence)}
      \item{\code{waterlog.50}}{Presence of waterlogging characteristics down to 50 cm (1: presence, 0: absence)}
      \item{\code{waterlog.100}}{Presence of waterlogging characteristics down to 100 cm (1: presence, 0: absence)}
      \item{\code{ph.0.10}}{Soil pH in 0-10 cm depth.}
      \item{\code{ph.10.30}}{Soil pH in 10-30 cm depth.}
      \item{\code{ph.30.50}}{Soil pH in 30-50 cm depth.}
      \item{\code{ph.50.100}}{Soil pH in 50-100 cm depth.}
      \item{\code{timeset}}{Factor with range of sampling year and label for sampling type for soil pH. no label: \eqn{CaCl_{2}}{CaCl_2} laboratory measurements, \code{field}: field estimate by indicator solution, \code{ptf}: \eqn{H_{2}0}{H_2_0} laboratory measurements transferred by pedotransfer function (univariate linear regression) to level of \eqn{CaCl_{2}}{CaCl_2} measures.}
       \item{\code{cl_mt_etap_pe}}{columns 14 to 238 contain environmental covariates representing soil forming factors. For more information see Details below.}
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
\strong{Soil data}

The soil data originates from various soil sampling campaigns since 1968. Most of the data was collected in conventional soil surveys in the 1970ties in the course of amelioration and farm land exchanges. As frequently observed in legacy soil data sampling site allocation followed a purposive sampling strategy identifying typical soils in an area in the course of polygon soil mapping.

\bold{\code{dclass}} contains drainage classes of three levels.
Swiss soil classification differentiates stagnic (I), gleyic (G) and anoxic/reduced (R) soil profile qualifiers with each 4, 6 resp. 5 levels. To reduce complexity the qualifiers  I, G and R  were aggregated to the degree of hydromorphic
	characteristic of a site with the ordered levels \emph{\code{well}} (qualifier labels I1--I2, G1--G3, R1 and no hydromorphic qualifier), \emph{\code{moderate}} well drained (I3--I4, G4) and \emph{\code{poor}} drained (G5--G6, R2--R5).


\bold{\code{waterlog}} indicates the \emph{\code{presence}} or \emph{\code{absence}} of waterlogging characteristics down 30, 50 and 100 cm soil depth. The responses were based on horizon qualifiers \sQuote{gg} or \sQuote{r} of the Swiss classification (\cite{Brunner et al. 1997}) as those were considered to limit plant growth. A horizon was given the qualifier \sQuote{gg} if it was strongly gleyic predominantly oxidized (rich in \eqn{Fe^{3+}}{Fe^3+}) and \sQuote{r} if it was anoxic predominantly reduced (poor in \eqn{Fe^{3+}}{Fe^3+}).


\bold{\code{pH}} was mostly sampled following genetic soil horizons. To ensure comparability between sites pH was transferred to fixed depth intervals of 0--10, 10--30, 30--50 and 50--100 cm by weighting soil horizons falling into a given interval. The data contains laboratory measurements that solved samples in \eqn{CaCl_{2}}{CaCl_2} or \eqn{H_{2}0}{H_2_0}. The latter were transferred to the level of \eqn{CaCl_{2}}{CaCl_2} measurements by univariate linear regression (label \code{ptf} in \emph{\code{timeset}}). Further, the dataset contains estimates of pH in the field by an indicator solution (Hellige pH, label \code{field} in \emph{\code{timeset}}).
The column \emph{\code{timeset}} can be used to partly correct for the long sampling period and the different sampling methods.


\strong{Environmental covariates}

The numerous covariates were assembled from the available spatial data in the case study area.
Each covariate name was given a prefix:

\itemize{
  \item \emph{\code{cl_}} climate covariates as precipitation, temperature, radiation
  \item \emph{\code{tr_}} terrain attributes, covariates derived from digital elevation models
  \item \emph{\code{ge_}} covariates from geological maps
  \item \emph{\code{sl_}} covariates from an overview soil map
}

References to the used datasets can be found in \cite{Nussbaum et al. 2017b}.

}
\references{
%%  ~~ possibly secondary sources and usages ~~
Brunner, J., Jaeggli, F., Nievergelt, J., and Peyer, K. (1997). Kartieren und Beurteilen von Landwirtschaftsboeden.
FAL Schriftenreihe 24, Eidgenoessische Forschungsanstalt fuer Agraroekologie und Landbau, Zuerich-Reckenholz (FAL).

Nussbaum, M., Spiess, K., Baltensweiler, A., Grob, U., Keller, A., Greiner, L., Schaepman, M. E., and Papritz, A., 2017b. Evaluation of digital soil mapping approaches with large sets of environmental covariates, SOIL Discuss., https://www.soil-discuss.net/soil-2017-14/, in review.

}
\examples{
data(berne)
}
\keyword{datasets}
