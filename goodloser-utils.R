## MJ
## 2017

## Theme for ggplot2.
theme_m <- function(...) {
  theme(text = element_text(size = 12,
                            colour = "black"),
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.line = element_line(),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(1, "mm"),
        panel.spacing.x = unit(7.5, "mm"),
        panel.spacing.y = unit(2, "mm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10,
                                  face = "bold"),
        strip.text.y = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom",
        ...)
}
## Mikael Poul Johannesson
## 2017

## Estimator that does cluster-robust standard errors.
estimator_regression <- function(formula, data,
                                 cluster = NULL,
                                 weights = NULL) {

  if (is.character(formula)) formula <- as.formula(formula)

  all_vars <- all.vars(formula)
  if (!is.null(cluster)) all_vars <- c(all_vars, cluster)
  if (!is.null(weights)) all_vars <- c(all_vars, weights)
  data <- data[, all_vars]
  data <- na.omit(data)

  if (!is.null(weights)) {
    we <- data[[weights]]
  } else {
    we <- NULL
  }

  fit <- lm(formula, data, weights = we)

  if (!is.null(cluster)) {

    cl <- data[[cluster]]
    M <- length(unique(cl))
    N <- length(cl)

    dfc <- (M / (M - 1)) * ((N - 1) / (N - fit$rank))
    u <- apply(estfun(fit), 2, function(x) tapply(x, cl, sum))
    vcov_cl <- dfc * sandwich(fit, meat = crossprod(u) / N)

    out <- coeftest(fit, vcov_cl)
    out <- tidy(out)

  } else {
    out <- tidy(fit)
  }

  return(out)
}

## Wrapper that returns AMCEs in a tidy manner.
amce <- function(data, post, treatments, subgroup = NULL, subset = NULL, cluster = NULL) {

  treat <- treatments
  for (var in c(treat, subgroup)) {
    if (!is.factor(data[[var]])) data[[var]] <- factor(data[[var]])
  }

  if (!is.null(subgroup)) {
    subgroup_values <- lapply(subgroup, function(x) levels(factor(data[[x]])))
    subgroup_values <- c(subgroup_values, list(treat))
    names(subgroup_values) <- c(subgroup, "treat")
    grid <- expand.grid(subgroup_values)
  } else {
    grid <- data.frame(treat = treat, stringsAsFactors = FALSE)
  }

  res <- lapply(1:nrow(grid), function(i) {

    component <- as.character(grid$treat[i])
    data_inner <- data

    if (!is.null(subgroup)) {
      for (j in 1:(ncol(grid) - 1)) {
        data_inner <- data_inner[as.character(data_inner[[names(grid)[j]]]) == as.character(grid[i, j]), ]
        if (nrow(data_inner) == 0) {
          warning(paste0("No rows in subgroup: ",
                         paste0(sapply(names(grid), function(x) paste0(x, "==", as.character(grid[[x]][i]))),
                                collapse = ", ")))
          return(NULL)
        }
      }
    }

    fit <- estimator_regression(formula = paste0(post, " ~ ", component),
                                data = data_inner,
                                cluster = cluster)

    baseline <- data.frame(term = levels(data[[component]])[1],
                           estimate = 0,
                           std.error = 0,
                           statistic = NA,
                           p.value = NA,
                           stringsAsFactors = FALSE)

    estimate <- fit[fit$term != "(Intercept)", ]
    estimate <- rbind(estimate, baseline)
    estimate$term <- gsub(paste0("^", component), "", estimate$term)
    estimate$term <- factor(estimate$term, levels = levels(factor(data_inner[[component]])))
    estimate <- estimate[order(estimate$term), ]
    estimate$term <- as.character(estimate$term)
    estimate$treatment <- component
    names(estimate)[names(estimate) == "term"] <- "value"

    if (!is.null(subgroup)) {
      for (j in 1:(ncol(grid) - 1)) {
        estimate[[names(grid)[j]]] <- as.character(grid[i, j])
      }
    }
    return(estimate)
  })

  res <- do.call("rbind", res)
  res$treatment <- factor(res$treatment, levels = unique(res$treatment))
  res$value <- factor(res$value, levels = unique(res$value))
  res$value_order <- 1:nrow(res)
  names(res) <- gsub("\\.", "_", names(res))
  res <- res[, c(subgroup, "treatment", "value", "value_order", "estimate", "std_error", "statistic", "p_value")]

  return(res)
}
