globalVariables("rowname")

#' Summary of confidence intervals object
#'
#' @param object Output from \code{\link{bootConfInt}}
#' @param ... Further arguments
#' @export
summary.BIGLconfInt <- function(object, ...) {

ans <- list()
ans$estimate = object$single$meanEffect
ans$sigLevel = paste0(round(object$cutoff*100), "%")
ans$singleCI = object$single$confIntMeanEffect
ans$call = object$single$Call

ans$confInt = object$offAxis[object$offAxis$call %in% c("Syn", "Ant"),]
ans$confInt[, c("estimate", "lower", "upper")] = round(ans$confInt[, c("estimate", "lower", "upper")], 4)
ans$totals <- data.frame("Syn" = sum(object$offAxis$call == "Syn"),
                         "Ant" = sum(object$offAxis$call == "Ant"),
                         "Total" = nrow(object$offAxis))
rownames(ans$totals) = ""

class(ans) <- append("summary.BIGLconfInt", class(ans))
ans
}

#' Print summary of BIGLconfInt object
#'
#' @param x Summary of BIGLconfInt object
#' @inheritParams summary.BIGLconfInt
#' @export
print.summary.BIGLconfInt <- function(x, ...) {

    #Overall
    cat("Overall effect\n")
    cat(sep = "", "Estimated mean departure from null response surface with ",
        x$sigLevel, " confidence interval:\n", round(x$estimate, 4), " [", round(x$singleCI[1], 4), ", ", round(x$singleCI[2], 4), "]\n")
    cat("Evidence for effects in data:", x$call, "\n\n")

    #Pointwise
    cat("Significant pointwise effects\n")
    print(x$confInt)
    cat("\nPointwise", x$sigLevel, "confidence intervals summary:\n")
    print(x$totals)
    cat("\n")
}

#' Plot confidence intervals in a contour plot
#'
#' @param x off axis confidence intervals, a data frame
#' @param color analysis with which to colour cells, either \code{effect-size} or \code{maxR}
#' @param showAll show all intervals in the plot or only significant ones, logical defaulting to \code{TRUE}
#' @param digits Numeric value indicating the number of digits used for numeric values
#' @param xlab String for the x axis label
#' @param ylab String for the y axis label
#' @param greyScale If \code{greyScale = TRUE}, then plot is in grey scale,
#'   otherwise in colour.
#' @param ... additional arguments, currently ignored
#' @importFrom stats setNames
#' @export
#' @note written after the contour() function in the \code{drugCombo} package
plot.BIGLconfInt <- function(x, color = "effect-size", showAll = TRUE, digits = 3, xlab, ylab, greyScale = FALSE, ...) {
  
  if (missing(xlab)) xlab <- sprintf("Dose (%s)", x$names[1])
  if (missing(ylab)) ylab <- sprintf("Dose (%s)", x$names[2])
  
  if ("maxR" %in% names(x)) {
    synOut <- x$maxR$Ymean
    names(synOut)[names(synOut) == "call"] <- "synCall"
    
    effectOut <- x$confInt$offAxis
    names(effectOut)[names(effectOut) == "call"] <- "effectCall"
    effectOut$d1 <- as.numeric(gsub("(.+)_.+", "\\1", rownames(effectOut)))
    effectOut$d2 <- as.numeric(gsub(".+_(.+)", "\\1", rownames(effectOut)))
    
    x <- merge(synOut, effectOut, by = c("d1","d2"))
  } else {
    x <- x$confInt$offAxis
    names(x)[names(x) == "call"] <- "effectCall"
    #show doses on equidistant grid
    d1d2 <- rownames(x)
    d1d2split <- sapply(d1d2, function(y) strsplit(y, split = "_")[[1]])
    x$d1 <- as.numeric(d1d2split[1,])
    x$d2 <- as.numeric(d1d2split[2,])
  }
  
  # prepare fill legend
  synCalls <- c("None", "Ant", "Syn")
  
  if (color == "effect-size") {
    x$synLabel <- factor(x$effectCall, labels = synCalls, levels = c("None", "Ant", "Syn"))
  } else {
    x$synLabel <- factor(x$synCall, labels = synCalls, levels = c("None", "Ant", "Syn"))
  }
  
  if(greyScale){
    legendColors <- c("grey70", "#636363", "#FEFCFF")
  } else {
    legendColors <- c("white", "pink", "lightblue")
  }
  
  names(legendColors) <- synCalls
  # subset to only the colors that are present in the data
  # legendColors <- legendColors[names(legendColors) %in% as.character(unique(x$synLabel))]
  
  # text to show
  fmt <- sprintf("%%.%if\n(%%.%if, %%.%if)", digits, digits, digits)
  if (isTRUE(showAll)) {
    x$label <- sprintf(fmt, x$estimate, x$lower, x$upper)
  } else {
    x$label <- ifelse(x$synLabel != "None",
                      sprintf(fmt, x$estimate, x$lower, x$upper),
                      "")
  }
  
  x$d1 <- factor(x$d1, levels = sort(unique(x$d1)),
                 labels = sort(unique(x$d1)), ordered = TRUE)
  x$d2 <- factor(x$d2, levels = sort(unique(x$d2)),
                 labels = sort(unique(x$d2)), ordered = TRUE)
  
  p <- ggplot(data = x, aes(x = .data$d1, y = .data$d2)) +
    geom_tile(aes(fill = .data$synLabel), color = "grey") +
    geom_text(aes(label = .data$label), show.legend = FALSE, size = 3) +
    # invisible points, used only for labels
    geom_point(aes(color = .data$synLabel), alpha = 0) +
    # round dose labels to digits
    scale_x_discrete(labels = format(as.numeric(levels(x$d1)), digits = digits)) +
    scale_y_discrete(labels = format(as.numeric(levels(x$d2)), digits = digits)) +
    scale_fill_manual(values = legendColors,
                      guide = "none", drop = FALSE) +
    scale_color_manual( # for a nicer legend
      values = setNames(1:3, nm = synCalls),
      limits = force,
      drop = FALSE,
      guide = guide_legend(title = "call:",
                           override.aes = list(alpha = 1, shape = 22, size = 8, color = "grey",
                                               fill = legendColors))
    ) +
    theme_minimal() +
    xlab(xlab) + ylab(ylab) +
    theme(
      panel.grid.major = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) 
  p
}


#' Plot confidence intervals from BIGL object in a contour plot
#'
#' @param BIGLobj Output from \code{\link{fitSurface}}
#' @param ... passed on to \code{\link{plot.BIGLconfInt}}
#' @export
plotConfInt <- function(BIGLobj, ...) {
  newBIGLobj <- BIGLobj
  class(newBIGLobj) <- ("BIGLconfInt")
  plot(newBIGLobj, ...)
}