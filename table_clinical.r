
tab.clin <- data.clinical_long %>% 
  ungroup() %>% 
  select(ID, numID, age, gender, Stroke.side,vol, clin.meas.name, clin.meas.value, tp) %>% 
  merge(data.GGP  %>% 
          ungroup() %>% 
          select(numID, tp, threshold, side, lab, GGP)
        ) %>% 
  filter(threshold == 1, ID %in% data.patients$ID & clin.meas.name %in% c('rGS','NIHSS', 'FM') & lab %in% c(1,3)) %>% 
  select(ID, age, gender, Stroke.side, vol, tp, clin.meas.name, clin.meas.value, lab, side, GGP) %>% 
  arrange(ID, tp)

tab.clin %>% 
  pivot_wider(id_cols = c(ID, tp)
              , names_from = c(clin.meas.name, lab)
              , values_from = c(clin.meas.value, GGP)
              ) %>% 
  flextable() %>% 
  merge_v()

tab.clin %>% 
  mutate(GGP.name=as_labeller(GGP_names)(lab) %>% as.character()) %>% 
  rename(GGP.value = GGP) %>%
  mutate(ID = as.numeric(factor(ID))) %>% gt
  dplyr::select(ID, age, gender, Stroke.side, vol, tp, clin.meas.name,clin.meas.value, side, GGP.name, GGP.value) %>% 
  write.csv(file = 'AppendixAdata.csv', sep = ',', row.names = FALSE, col.names = TRUE, na = 'NA')

tab.clin %>% 
  gather(variable, value, GGP) %>% 
  unite(temp, lab, side, variable) %>% 
  spread(temp, value) %>%
  gather(variable, value, clin.meas.value) %>% 
  unite(temp, clin.meas.name, variable) %>% 
  spread(temp, value) %>%
  select(ID, age, gender, Stroke.side, vol, tp
         , NIHSS_clin.meas.value, FM_clin.meas.value, rGS_clin.meas.value
         , `1_Left_GGP`, `3_Left_GGP`, `1_Right_GGP`, `3_Right_GGP`) %>% 
  flextable() %>% 
  add_header_row(values = c('ID', 'Age', 'Sex','Stroke side', 'Lesion volume', 'Time'
                            , 'NIHSS', 'FM', 'rGS'
                            , 'Left Hemisphere', 'Right Hemisphere'
                            )
                 , colwidths = c(1,1,1,1,1,1,1,1,1,2,2)) %>% 
  set_header_labels(values = list(ID = 'ID'
                                  , age = 'Age'
                                  , gender = 'Sex'
                                  , Stroke.side = 'Stroke side'
                                  , vol = 'Lesion volume'
                                  , tp = 'Time'
                                  , `1_Left_GGP` = 'Efficiency'
                                  , `1_Right_GGP` = 'Efficiency'
                                  , `3_Left_GGP` = 'Modularity'
                                  , `3_Right_GGP` = 'Modularity'
                                  , `FM_clin.meas.value` = 'FM'
                                  , `NIHSS_clin.meas.value` = 'NIHSS'
                                  , `rGS_clin.meas.value` = 'rGS'
                                  )) %>% 
  merge_v(part = 'body', j = c(1,2,5,6,10,11,12)) %>% 
  merge_v(part = 'header') %>%
  flextable::theme_vanilla()
