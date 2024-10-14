

library(tidymodels)
library(skimr)
library(iml)
library(rpart.plot)
library(fastshap)
library(ggbeeswarm)

winequality <- readr::read_csv()#读取数据
colnames(winequality)

for(i in c(12)){
  winequality[[i]] <- factor(winequality[[i]])
}#将第i个变量转换为factor(分类变量):c(i)

skimr::skim(winequality)#数据概况

set.seed(4321)
datasplit <- initial_split(winequality, prop = 0.75, strata = quality)#数据拆分
traindata <- training(datasplit)
testdata <- testing(datasplit)

#数据预处理定义：'quality ~ .' 响应变量为quality，其他全为预测变量
datarecipe <- recipe(quality ~ ., traindata) %>%
  step_naomit(all_predictors()) %>%  #删除预测变量中的NA值
  prep()
#datarecipe

traindata2 <- bake(datarecipe, new_data = NULL) %>%
  dplyr::select(quality, everything())
testdata2 <- bake(datarecipe, new_data = testdata) %>%
  dplyr::select(quality, everything())

skimr::skim(traindata2)
skimr::skim(testdata2)




#模型定义
model_dt <- decision_tree( #决策树
  mode = "classification", #分类
  engine = "rpart",        #决策树引擎
  tree_depth = tune(),     #最大深度
  min_n = tune(),          #每个节点的最小样本数
  cost_complexity = tune() #树的复杂度
) %>%
  set_args(model=TRUE)     #模型额外参数（工作流中使用全部参数）

model_dt

# workflow (tidymodels工作流)
wk_dt <- 
  workflow() %>%
  add_model(model_dt) %>%
  add_formula(quality ~ .)
wk_dt

# 重抽样设定-5折交叉验证
set.seed(42)
folds <- vfold_cv(traindata2, v = 5)
folds

# 超参数寻优范围
set.seed(42)
hpset_dt <- parameters(tree_depth(range = c(3, 7)),min_n(range = c(5, 10)),cost_complexity(range = c(-5, -1)))#定义超参数设置
#hpgrid_dt <- grid_regular(hpset_dt, levels = c(3, 2, 3))#创建规则参数网格
hpgrid_dt <- grid_random(hpset_dt, size = 10)#创建随机参数网格
hpgrid_dt
log10(hpgrid_dt$cost_complexity)


# 交叉验证网格搜索过程
set.seed(42)
tune_dt <- wk_dt %>%
  tune_grid(resamples = folds,
            grid = hpgrid_dt,
            metrics = metric_set(yardstick::accuracy, 
                                 yardstick::roc_auc,
                                 yardstick::pr_auc),
            control = control_grid(save_pred = T, verbose = T))

# 图示交叉验证结果
autoplot(tune_dt)
eval_tune_dt <- tune_dt %>%
  collect_metrics()
eval_tune_dt

# 经过交叉验证得到的最优超参数
hpbest_dt <- tune_dt %>%
  select_by_one_std_err(metric = "accuracy", desc(cost_complexity))
hpbest_dt

# 采用最优超参数组合训练最终模型
final_dt <- wk_dt %>%
  finalize_workflow(hpbest_dt) %>%
  fit(traindata2)
final_dt

# 提取最终的算法模型
final_dt2 <- final_dt %>%
  extract_fit_engine()

# 树形图
rpart.plot(final_dt2)

# 变量重要性
final_dt2$variable.importance
# par(mar = c(10, 5, 1, 1))
barplot(final_dt2$variable.importance, las = 2)

#################################################################

# 应用模型-预测训练集
predtrain_dt <- final_dt %>%
  predict(new_data = traindata2, type = "prob") %>%
  bind_cols(final_dt %>%
              predict(new_data = traindata2, type = "class")) %>%
  bind_cols(traindata2 %>% select(quality)) %>%
  mutate(dataset = "train")
# 将预测分类的levels设定成因变量实际值的levels
levels(predtrain_dt$.pred_class) <- levels(predtrain_dt$quality)
predtrain_dt

# 评估模型ROC曲线-训练集上
roctrain_dt <- predtrain_dt %>%
  roc_curve(quality, .pred_H:.pred_M) %>%
  mutate(dataset = "train")
roctrain_dt
autoplot(roctrain_dt)

# 混淆矩阵
cmtrain_dt <- predtrain_dt %>%
  conf_mat(truth = quality, estimate = .pred_class)
cmtrain_dt
autoplot(cmtrain_dt, type = "heatmap") +
  scale_fill_gradient(low = "white", high = "skyblue") +
  theme(text = element_text(size = 15))
# 合并指标
eval_train_dt <- cmtrain_dt %>%
  summary() %>%
  bind_rows(predtrain_dt %>%
              roc_auc(quality, .pred_H:.pred_M)) %>%
  mutate(dataset = "train")
eval_train_dt

#################################################################结果不行返回调参

# 应用模型-预测测试集
predtest_dt <- final_dt %>%
  predict(new_data = testdata2, type = "prob") %>%
  bind_cols(final_dt %>%
              predict(new_data = testdata2, type = "class")) %>%
  bind_cols(testdata2 %>% select(quality)) %>%
  mutate(dataset = "test") %>%
  mutate(model = "dt")
# 将预测分类的levels设定成因变量实际值的levels
levels(predtest_dt$.pred_class) <- levels(predtest_dt$quality)
predtest_dt
# 评估模型ROC曲线-测试集上
roctest_dt <- predtest_dt %>%
  roc_curve(quality, .pred_H:.pred_M) %>%
  mutate(dataset = "test")
autoplot(roctest_dt)


# 混淆矩阵
cmtest_dt <- predtest_dt %>%
  conf_mat(truth = quality, estimate = .pred_class)
cmtest_dt
autoplot(cmtest_dt, type = "heatmap") +
  scale_fill_gradient(low = "white", high = "skyblue") +
  theme(text = element_text(size = 15))
# 合并指标
eval_test_dt <- cmtest_dt %>%
  summary() %>%
  bind_rows(predtest_dt %>%
              roc_auc(quality, .pred_H:.pred_M)) %>%
  mutate(dataset = "test")
eval_test_dt

#################################################################

# 合并训练集和测试集上ROC曲线
roctrain_dt %>%
  bind_rows(roctest_dt) %>%
  mutate(dataset = factor(dataset, levels = c("train", "test"))) %>%
  ggplot(aes(x = 1-specificity, y = sensitivity, color = dataset)) +
  geom_path(linewidth = 1) +
  facet_wrap(~.level) +
  theme_bw()

# 合并训练集和测试集上性能指标
eval_dt <- eval_train_dt %>%
  bind_rows(eval_test_dt) %>%
  mutate(model = "dt")
eval_dt

#############################################################

# 最优超参数的交叉验证指标平均结果
eval_best_cv_dt <- eval_tune_dt %>%
  inner_join(hpbest_dt[, 1:3])
eval_best_cv_dt

# 最优超参数的交叉验证指标具体结果
eval_best_cv5_dt <- tune_dt %>%
  collect_predictions() %>%
  inner_join(hpbest_dt[, 1:3]) %>%
  group_by(id) %>%
  roc_auc(quality, .pred_H:.pred_M) %>%
  ungroup() %>%
  mutate(model = "dt") %>%
  inner_join(eval_best_cv_dt[c(4,6,8)])
eval_best_cv5_dt


# 最优超参数的交叉验证指标图示
eval_best_cv5_dt %>%
  filter(.metric == "roc_auc") %>%
  ggplot(aes(x = id, y = .estimate, group = 1)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "", y = "roc_auc") +
  theme_bw()


# 最优超参数的交叉验证图示
tune_dt %>%
  collect_predictions() %>%
  inner_join(hpbest_dt[, 1:3]) %>%
  group_by(id) %>%
  roc_curve(quality, .pred_H:.pred_M) %>%
  ungroup() %>%
  ggplot(aes(x = 1-specificity, y = sensitivity, color = id)) +
  geom_path(linewidth = 1) +
  facet_wrap(~.level) +
  theme_bw()

###################################################################

# 自变量数据集
colnames(traindata2)
traindatax <- traindata2[,-1]#选择用于预测的变量
colnames(traindatax)



# 变量重要性-基于置换
predictor_model <- Predictor$new(
  final_dt, 
  data = traindatax,
  y = traindata2$quality
)
imp_model <- FeatureImp$new(
  predictor_model, 
  loss = "ce"
)
# 数值
imp_model$results
# 图示
imp_model$plot() +
  theme_bw()


# 变量效应
predictor_model <- Predictor$new(
  final_dt, 
  data = traindatax,
  y = traindata2$quality,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>% pull(1)
  }
)
pdp_model <- FeatureEffect$new(
  predictor_model, 
  feature = "pH",
  method = "pdp"
)
# 数值
pdp_model$results
# 图示
pdp_model$plot() +
  theme_bw()

# 所有变量的效应全部输出
effs_model <- FeatureEffects$new(predictor_model, method = "pdp")
# 数值
effs_model$results
# 图示
effs_model$plot()

# 单样本shap分析
shap_model <- Shapley$new(
  predictor_model, 
  x.interest = traindatax[1,]
)
# 数值
shap_model$results
# 图示
shap_model$plot() +
  theme_bw()

# 基于所有样本的shap分析
# fastshap包
shap <- explain(
  final_dt, 
  X = as.data.frame(traindatax),
  nsim = 10,
  adjust = T,
  pred_wrapper = function(model, newdata) {
    predict(model, newdata, type = "prob") %>% pull(1)
  }
)


data1 <- shap %>%
  as.data.frame() %>%
  dplyr::mutate(id = 1:n()) %>%
  pivot_longer(cols = -(ncol(traindatax)+1), values_to = "shap")
shapimp <- data1 %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(shap.abs.mean = mean(abs(shap))) %>%
  dplyr::arrange(shap.abs.mean) %>%
  dplyr::mutate(name = forcats::as_factor(name))
data2 <- traindatax  %>%
  dplyr::mutate(id = 1:n()) %>%
  pivot_longer(cols = -(ncol(traindatax)+1))

# 所有变量shap图示
data1 %>%
  left_join(data2) %>%
  dplyr::rename("feature" = "name") %>%
  dplyr::group_by(feature) %>%
  dplyr::mutate(
    value = (value - min(value)) / (max(value) - min(value)),
    feature = factor(feature, levels = levels(shapimp$name))
  ) %>%
  dplyr::arrange(value) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = shap, y = feature, color = value)) +
  geom_quasirandom(width = 0.2) +
  scale_color_gradient(
    low = "red", 
    high = "blue", 
    breaks = c(0, 1), 
    labels = c(" Low", "High "), 
    guide = guide_colorbar(barwidth = 1, 
                           barheight = 20,
                           ticks = F,
                           title.position = "right",
                           title.hjust = 0.5)
  ) +
  labs(x = "SHAP value", color = "Feature value") +
  theme_bw() +
  theme(legend.title = element_text(angle = -90))


# 单变量shap图示
data1 %>%
  left_join(data2) %>%
  dplyr::rename("feature" = "name") %>%
  dplyr::filter(feature == "alcohol") %>%
  ggplot(aes(x = value, y = shap)) +
  geom_point() +
  geom_smooth(se = F, span = 0.5) +
  labs(x = "alcohol") +
  theme_bw()






