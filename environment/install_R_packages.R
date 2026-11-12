# Install the exact package versions used for the published analysis (R 4.4.2).
install.packages(c('data.table','remotes'))
remotes::install_github('MRCIEU/TwoSampleMR@6268c2cf74da8533a89094dda78c24db8e87bb25')   # 0.6.15
remotes::install_github('rondolab/MR-PRESSO@3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87')   # 1.0
remotes::install_github('qingyuanzhao/mr.raps@27b96f31e26ae97a395422bd757f514a7f96cc85') # 0.4.2
