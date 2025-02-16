# fully internal function
# given a data.frame try to infer
# what kind of meta-analysis data we're dealing with

detect_input_type <- function(data,
                              group="group",
                              treatment="treatment",
                              outcome="outcome") {
  # if(class == "baggr_data")
  if(!("data.frame" %in% class(data)))
    stop("Can't detect input type because it's not data.frame")

  # Summary data -----
  if("tau" %in% names(data) && "se" %in% names(data)){
    if("treatment" %in% names(data))
      return("pool_narrow") #will need conversion
    else
      return("pool_noctrl_narrow")
  }

  if(!any(is.na(match(c("tau", "mu", "se.mu", "se.tau"),
                      names(data)))))
    return("pool_wide")

  if (!any(is.na(match(c("interval_left", "interval_right"), 
                       names(data))))) {
    return("individual_censored")
                       }

  if (length(outcome) == 2) {
    return("individual_censored")
  }

  if(check_columns_binary(data, stop = FALSE))
    return("pool_binary")

  # Individual-level data -----
  if(!is.null(data[[group]])){
    if(nrow(data) > length(unique(data[[group]]))){
      if(is_binary(data[[outcome[1]]]))
        return("individual_binary")
      else
        return("individual")
    }
  }

  # Shrug emoji -----
  return("unknown")
}
