

som

class(sce)


features <- .get_features(x, features)
if (is.null(marker_classes(x))) {
  rowData(x)$marker_class <- factor(c("state", "type")[as.numeric(rownames(x) %in% 
                                                                    features) + 1], levels = c("type", "state", "none"))
}

fsom <- FlowSOM::ReadInput(flowFrame(t(assay(sce, "exprs"))))
FlowSOM::NewData(fsom, toTransform=markers_to_use)
