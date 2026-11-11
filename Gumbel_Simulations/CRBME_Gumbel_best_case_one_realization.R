# =====================================================================
# CR-BME GUMBEL FIGURES - ONE REALIZATION
# Selected case: l=3, q=(0.75,0.95)
# Same simulation model and tuning ranges as the original Gumbel study.
# A dedicated JOINT realization at 625 sites + 900 grid nodes is used.
# This is NOT one of the original 5 study realizations; their averages
# should not be reported as results from this figure.
# =====================================================================
library(devtools)
library(gstat)
library(BMEmapping)
library(sf)
library(ggplot2)
library(dplyr)
library(caret)

load_all()
dir.create('figs', showWarnings=FALSE)

# 1. Original spatial design and covariance model
set.seed(1234)
n <- 625
coords <- data.frame(x=runif(n,0,100), y=runif(n,0,100))
nh <- 1:125; ns <- 126:525; nk <- 526:625
vgm_model <- vgm(psill=0.5, model='Exp', range=60, nugget=0.5)
g.dummy <- gstat(formula=z~1, dummy=TRUE, beta=0,
                 model=vgm_model, nmax=50)
gum_mu <- 0.5
gum_beta <- 1
to_gumbel <- function(z) {
  p <- pmin(pmax(pnorm(z), .Machine$double.eps), 1-.Machine$double.eps)
  gum_mu - gum_beta*log(-log(p))
}
dgumbel_custom <- function(z) {
  t <- (z-gum_mu)/gum_beta
  exp(-t-exp(-t))/gum_beta
}

# 2. Chosen case and same 10 PAIRED CV settings as original study
l_fig <- 3
q_fig <- c(0.75,0.95)
n_grid <- 30
param_grid <- data.frame(lambda=seq(0.5,5,length.out=10),
                         mu=seq(1,5,length.out=10))

# 3. Jointly simulate at data sites and grid nodes (one realization).
# Separate seed intentionally produces a dedicated illustrative field.
grid_xy <- expand.grid(x=seq(0,100,length.out=n_grid),
                       y=seq(0,100,length.out=n_grid))
both <- rbind(coords,grid_xy)
both_sf <- st_as_sf(both,coords=c('x','y'))
set.seed(2026)
sim_both <- predict(g.dummy,newdata=both_sf,nsim=1)
z_field <- to_gumbel(sim_both$sim1)
z_data <- z_field[seq_len(n)]
z_truth <- z_field[n+seq_len(nrow(grid_xy))]

ch <- coords[nh,,drop=FALSE]; zh <- z_data[nh]
cs <- coords[ns,,drop=FALSE]; zs <- z_data[ns]
ck <- coords[nk,,drop=FALSE]; zk <- z_data[nk]

# 4. Soft intervals (same formula and width/placement settings)
ln <- runif(length(zs),min=l_fig-0.5,max=l_fig+0.5)
u  <- runif(length(zs),min=q_fig[1],max=q_fig[2])
a <- zs-ln*u
b <- zs+ln*(1-u)
data_object <- bme_map(ch,cs,zh,a,b)
zk_rng <- base::range(zh,a,b)

# 5. Fit CBME variogram using the original fitting procedure.
# Inspect fit and warnings; positional extraction assumes two-row fit.
df_vg <- st_as_sf(data.frame(ch,z=zh),coords=c('x','y'))
vg <- variogram(z~1,data=df_vg)
vgf <- fit.variogram(vg,model=vgm(c('Exp','Sph')))
print(vgf)
if (nrow(vgf)<2 || !all(is.finite(vgf$psill)) ||
    !all(is.finite(vgf$range))) {
  stop('Check fitted variogram before predicting.')
}
model_type <- as.character(vgf[2,1])
nugget <- vgf[1,2]
sill <- vgf[2,2]
range_par <- vgf[2,3]

# 6. Five-fold CV for CR-BME, tuning posterior MEAN with RMSE
set.seed(1234)
folds <- createFolds(zh,k=5)
cv_summary <- vector('list',nrow(param_grid))
for(g in seq_len(nrow(param_grid))) {
  fold_results <- vector('list',length(folds))
  for(i in seq_along(folds)) {
    test_idx <- folds[[i]]
    data_cv <- bme_map(ch=ch[-test_idx,,drop=FALSE], cs=cs,
                       zh=zh[-test_idx], a=a,b=b)
    preds <- new_bme_predict(x=ch[test_idx,,drop=FALSE], data=data_cv,
                   lambda=param_grid$lambda[g],mu=param_grid$mu[g],type='mean')
    if(length(preds$mean)!=length(test_idx) ||
       any(!is.finite(preds$mean))) stop('Non-finite CV prediction.')
    fold_results[[i]] <- data.frame(
      observed=zh[test_idx], predicted=preds$mean)
  }
  combined <- do.call(rbind,fold_results)
  residual <- combined$observed-combined$predicted
  cv_summary[[g]] <- data.frame(
    lambda=param_grid$lambda[g], mu=param_grid$mu[g],
    MAE=mean(abs(residual)),RMSE=sqrt(mean(residual^2)))
}
cv_table <- do.call(rbind,cv_summary)
best_row <- cv_table[which.min(cv_table$RMSE),,drop=FALSE]
best_lambda <- best_row$lambda
best_mu <- best_row$mu
print(best_row)
write.csv(cv_table,'figs/gumbel_cv_tuning.csv',row.names=FALSE)

# 7. Predict at grid nodes and 100 independent hold-out sites
predict_cbme <- function(x) bme_estimate(
  x,data_object,model_type,nugget,sill,range_par,zk_range=zk_rng)[,2]
predict_crbme <- function(x) new_bme_estimate(
  x,data_object,zk_range=zk_rng,
  lambda=best_lambda,mu=best_mu)[,2]
message('Predicting over ',nrow(grid_xy),' grid nodes')
pred_cbme <- predict_cbme(grid_xy)
pred_crbme <- predict_crbme(grid_xy)
pred_cbme_test <- predict_cbme(ck)
pred_crbme_test <- predict_crbme(ck)
if(any(!is.finite(c(pred_cbme,pred_crbme,pred_cbme_test,pred_crbme_test)))) {
  stop('Non-finite predictions. Review fitted variogram and model.')
}

# 8. Common theme
th <- theme_minimal(base_size=13)+theme(
  panel.border=element_rect(colour='grey30',fill=NA),
  panel.grid.minor=element_blank(),
  strip.text=element_text(face='bold',size=12),legend.position='right')

# Figure 1 - locations
loc_df <- coords %>% mutate(type=factor(
  c(rep('Hard (125)',length(nh)), rep('Soft (400)',length(ns)),
    rep('Hold-out (100)',length(nk))),
  levels=c('Hard (125)','Soft (400)','Hold-out (100)')))
p1 <- ggplot(loc_df,aes(x,y,shape=type,colour=type))+
  geom_point(size=2.1,stroke=0.8)+
  scale_shape_manual(values=c(16,1,4))+
  scale_colour_manual(values=c('#004886','#8AA7C4','#DC7828'))+
  coord_fixed(xlim=c(0,100),ylim=c(0,100))+
  labs(x='x',y='y',shape=NULL,colour=NULL)+th
ggsave('figs/gumbel_locations.pdf',p1,width=6.4,height=5)

# Figure 2 - empirical density of this realization vs true Gumbel pdf
zz <- seq(min(z_data),max(z_data),length.out=400)
ref <- data.frame(z=zz,d=dgumbel_custom(zz))
p2 <- ggplot()+
  geom_histogram(data=data.frame(z=z_data),
    aes(z,after_stat(density)),bins=40,fill='grey85',colour='white')+
  geom_density(data=data.frame(z=z_data),aes(z),
    colour='#8AA7C4',linewidth=0.8)+
  geom_line(data=ref,aes(z,d),colour='#DC7828',linewidth=1.2)+
  labs(title='Gumbel marginal density',subtitle='One joint spatial realization',
       x='z',y='Density')+th
ggsave('figs/gumbel_density.pdf',p2,width=6.4,height=4.4)

# Figure 3 - actual simulated truth at every grid node and predictions
map_df <- rbind(
  data.frame(grid_xy,value=z_truth,panel='Truth'),
  data.frame(grid_xy,value=pred_cbme,panel='CBME'),
  data.frame(grid_xy,value=pred_crbme,panel='CR-BME'))
map_df$panel <- factor(map_df$panel,levels=c('Truth','CBME','CR-BME'))
p3 <- ggplot(map_df,aes(x,y,fill=value))+
  geom_raster(interpolate=FALSE)+facet_wrap(~panel,nrow=1)+
  scale_fill_viridis_c(option='D',limits=range(map_df$value),name='z')+
  coord_fixed(expand=FALSE)+
  labs(title='Gumbel spatial prediction',
       subtitle='One realization; l = 3, q = (0.75, 0.95)',x='x',y='y')+th
ggsave('figs/gumbel_prediction_maps.pdf',p3,width=11,height=4.5)

# Figure 4 - error maps from exactly simulated truth at grid nodes
err_df <- rbind(
  data.frame(grid_xy,err=abs(pred_cbme-z_truth),panel='CBME'),
  data.frame(grid_xy,err=abs(pred_crbme-z_truth),panel='CR-BME'))
p4 <- ggplot(err_df,aes(x,y,fill=err))+
  geom_raster(interpolate=FALSE)+facet_wrap(~panel,nrow=1)+
  scale_fill_viridis_c(option='magma',limits=range(err_df$err),
                       name='|error|')+
  coord_fixed(expand=FALSE)+
  labs(title='Gumbel absolute grid prediction errors',
       subtitle='One realization; l = 3, q = (0.75, 0.95)',x='x',y='y')+th
ggsave('figs/gumbel_error_maps.pdf',p4,width=8.5,height=4.5)

# 9. Display accuracy on grid AND at original 100 holdout sites
metric_row <- function(method,truth,pred,site_type) data.frame(
  Evaluation=site_type,Method=method,
  MAE=mean(abs(truth-pred)),MSE=mean((truth-pred)^2),
  RMSE=sqrt(mean((truth-pred)^2)))
results <- rbind(
  metric_row('CBME',z_truth,pred_cbme,'900 grid nodes'),
  metric_row('CR-BME',z_truth,pred_crbme,'900 grid nodes'),
  metric_row('CBME',zk,pred_cbme_test,'100 hold-out sites'),
  metric_row('CR-BME',zk,pred_crbme_test,'100 hold-out sites'))
print(results)
write.csv(results,'figs/gumbel_accuracy.csv',row.names=FALSE)
cat('\nCV-selected lambda:',best_lambda,'mu:',best_mu,'\n')
print(p1); print(p2); print(p3); print(p4)
