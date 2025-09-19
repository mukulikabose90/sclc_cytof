library(ROCR)
library(caret)
source("source/cytof_de_function.R")

set.seed(42)
################################################################################

sce <- readRDS("data/cytof_objects/sclc_cytof_sce_object_with_clusters_3.rds")

sce <- sce[rowData(sce)$marker_class == "state",]

y <- assay(sce,"exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

df$result <- ifelse(df$new_clusters %in% c(1,6,10), 1, 0)

data <- df[,c(1:16,31)]


smp_size <- floor(0.75 * nrow(data))

train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]

model <- glm("result ~ .", data=train, family = "binomial")

summary(model)


logitModelPred <- predict(model, test, type = "response")

# setting the cut-off probablity
classify50 <- ifelse(logitModelPred > 0.5,"1","0")

# ordering the levels
classify50 <- ordered(classify50, levels = c("1", "0"))
test$result <- ordered(test$result, levels = c("1", "0"))

# confusion matrix
cm <- table(Predicted = classify50, Actual = test$result)
cm

###############
data$result <- as.factor(data$result)


Train <- createDataPartition(data$result, p=0.6, list=FALSE)
training <- data[ Train, ]
testing <- data[ -Train, ]

ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE)



mod_fit <- train(result ~ CD24 + EpCAM + ASCL1 + Twist + 
                   POU2F3 +`p-YAP`+ DLL3 + Vimentin + NeuroD1,  data=training, method="glm", family="binomial",
                 trControl = ctrl, tuneLength = 5)

mod_fit <- train(result ~ .,  data=training, method="glm", family="binomial",
                 trControl = ctrl, tuneLength = 5)

pred = predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$result)

######################
confusionMatrix(cm)


PredLR <- predict(model, test,type = "response")
lgPredObj <- prediction((1-PredLR),test$result)
lgPerfObj <- performance(lgPredObj, "tpr","fpr")
# plotting ROC curve
plot(lgPerfObj,main = "ROC Curve",col = 2,lwd = 2)
abline(a = 0,b = 1,lwd = 2,lty = 3,col = "black")



# area under curve
aucLR <- performance(lgPredObj, measure = "auc")
aucLR <- aucLR@y.values[[1]]
aucLR

