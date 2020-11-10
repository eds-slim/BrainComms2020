d <- data.clinical_long %>%
  mutate(clin.meas.value = if_else(clin.meas.name == "FM", 66 - clin.meas.value, clin.meas.value)) %>%
  merge(d.LGP) %>%
  filter(clin.meas.name %in% c("NIHSS", "FM", "rGS") & relpos == "ipsi" & tp >= 1 & numID != 40 &
    LGP.name %in% c("strength", "efficiency", "clustering"))


d.LGP.clinical.joint <- d %>%
  group_by(node, LGP.name, clin.meas.name) %>%
  nest() %>%
  # mutate(plt= pmap(list(data,LGP.name,clin.meas.name)
  #                 ,~joint_fit_CV('LGP.value.diff',..1,..2,..3
  #                                , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
  #                 )
  # )
  # ) %>%
  mutate(mdl = pmap(
    list(data, LGP.name, clin.meas.name),
    ~ joint_fit_CV_mdl("LGP.value.diff", ..1, ..2, ..3,
      fam = if_else(..3 %in% c("NIHSS", "FM"), "quasipoisson", "gaussian")
    )
  ))



# p.LGP.clinical.joint <- do.call(arrangeGrob,c(d.LGP.clinical.joint$plt
#                                              , ncol=length(unique(d.LGP.clinical.joint$clin.meas.name))
#                                              , as.table=F)
# )




tbl.LGP.clinical <- d.LGP.clinical.joint %>%
  mutate(tidy = map(mdl, tidy)) %>%
  unnest(tidy) %>%
  filter(term == "fitted(fm)") %>%
  group_by(node, LGP.name, clin.meas.name) %>%
  nest() %>%
  mutate(p.min = map_dbl(data, ~ min(.$p.value))) %>%
  filter(p.min < 0.05 / 36) %>%
  dplyr::select(-p.min) %>%
  unnest(data) %>%
  mutate(stats = cell_spec(sprintf("%1.2f +/- %1.2f, p=%1.4f", estimate, std.error, p.value)
                           , format = "html"
                           , bold = p.value < 0.05)) %>%
  dplyr::select(-c(term, estimate, std.error, statistic, p.value)) %>%
  unite(temp, clin.meas.name, LGP.name) %>%
  spread(temp, stats, fill = "") %>%
  kable(format = "html", escape = F) %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  collapse_rows(columns = 1)

tbl.LGP.clinical
