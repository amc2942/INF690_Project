library(raster)
library(ggplot2)

statusquo_stack = rast('statusquo/mso_suitability_raster_stack.tif')
fast_4fri_stack = rast('fast4fri/mso_suitability_raster_stack.tif')

# Make up a "Good habitat" threshold using 75th percentile of habitat quality in year 1
year1_suitability_values = values(statusquo_stack[[1]])
thresh = quantile(year1_suitability_values[year1_suitability_values>0],na.rm = T,probs=c(.75))[[1]]
# Add up area greater than threshold in each year
mso_suitable_area_fast = apply(values(fast_4fri_stack>=thresh),FUN=sum,MARGIN=2)
mso_suitable_area_statusquo = apply(values(statusquo_stack>=thresh),FUN=sum,MARGIN=2)
mso_suitable_area = data.table(Year = seq(1:length(mso_suitable_area_fast)),
                               MSOSuitableAreaFast4FRI=mso_suitable_area_fast,
                               MSOSuitableAreaStatusQuo=mso_suitable_area_statusquo
)

ggplot(mso_suitable_area[1:cutoff_year]) +
  geom_line(aes(Year,MSOSuitableAreaFast4FRI,color='Fast 4FRI')) +
  geom_line(aes(Year,MSOSuitableAreaStatusQuo,color='Status Quo'))+
  labs(title = 'Suitable Habitat for Mexican Spotted Owl',color='Scenario')+
  ylab('Hectares of suitable habitat') +
  xlab('Year of Simulation')

ggsave('mso_suitability_trend_comparison.png')