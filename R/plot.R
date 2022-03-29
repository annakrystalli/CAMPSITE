#' Plot CAMPSITE simulation results
#'
#' @param x results object of class "cs_sim_results" output from `sim_CAMPSITE`.
#' @param ylims numeric vector of length 2 specifying limits for the y axis. If 
#' `NULL` (default), limits are determined automatically
#'
#' @return plot of simulated lineage evolution across trait space
#' @export
#'
#' @examples
#' \dontrun{
#' res <- cs_simulate(pars = cs_pars(), ou = list(opt = NULL, alpha4 = NULL), root.value = 0, age.max = 10, 
#'                    age.ext = NULL, step_size = 0.01, bounds = c(-Inf, Inf), 
#'                    plot = FALSE, ylims = NULL, full_results = TRUE) 
#' 
#' plot(res)                   
#' }
#'
plot.cs_sim_results <- function(x, ylims = NULL) {
  if(!x$full_results){
    stop("No lineage and trait data available to plot")
  }
  cs_plot_trait_evolution(x$traits, x$lineages, x$t_end, x$step_size, 
                          ylims = ylims)
}


cs_plot_trait_evolution <- function(traits, lineages, t, step_size, ylims = NULL) {
  # plots a simulation, with incipient lineages in red, good in black
  lineages[lineages[, "end_time"] < 0, "end_time"] <- t
  
  steps <- t / step_size
  # extract maximal trait values for each lineage
  max_trait <- sapply(traits, max_trait_value) 
  # extract minimal trait values for each lineage
  min_trait <- sapply(traits, min_trait_value) 
  
  
  if (is.null(ylims)) {
    # auto generate plot area
    plot(1,
         type = "n",
         xlim = c(1, steps),
         ylim = c(min(min_trait), 
                  max(max_trait)),
         ylab = "trait value",
         xlab = paste0("Time steps (step size = ", step_size ,")")
    ) 
  } else {
    # generate plot area with pre-determined limits for trait space
    plot(1, type = "n",
         xlim = c(1, steps),
         ylim = ylims, 
         ylab = "trait value",
         xlab = paste0("Time steps (step size = ", step_size ,")")
    )
    
  }
  
  opacity <- 0.6
  colg <- scales::alpha("black", opacity)
  coli <- scales::alpha("red", opacity)
  
  
  completion <- lineages[, "spec_or_ext_ct"] #extract vector of speciation or extinction times for each lineage
  #handle and plot each lineage based on its end point:
  for (i in seq_along(completion)) {
    
    trait_values <- extract_trait_values(traits[[i]])
    
    if (is.na(completion[i])) { #extinct incipient lineages
      graphics::lines(x = seq_along(trait_values), # branch length
            y = trait_values, #trait values
            col = coli, #incipient lineages in red
            cex = 0.5)
    }
    else if (lineages[i, "status"] == -2 &&
             !is.na(lineages[i, "spec_ct"]) &&
             completion[i] != lineages[i, "spec_ct"]) { # extinct good lineages
      spec_t_s <- lineages[i, "spec_ct"] / step_size + 1
      
      graphics::lines(x = seq(1, spec_t_s), # branch length of incipient stage
            y = trait_values[seq(1, spec_t_s)],
            col = coli, # incipient stage in red
            cex = 0.5)
      graphics::lines(x = seq(spec_t_s, length(trait_values)), #branch length of good stage
            y =  trait_values[seq(spec_t_s, length(trait_values))],
            col = colg, # good stage in black
            cex = 0.5)
    }
    else { #  extant good and incipient lineages
      comp_t_s <- completion[i] / step_size + 1
      graphics::lines(x = seq(1, comp_t_s), #branch length of incipient stage
            y = trait_values[seq(1, comp_t_s)], #trait values of incipient stage
            col = coli, #incipient stage in red
            cex = 0.5)
      if (length(trait_values) > round(comp_t_s)) { #if completed speciation before end of process
        graphics::lines(x = seq(comp_t_s, length(trait_values)), #branch length of good stage
              y = trait_values[seq(comp_t_s, length(trait_values))], #trait values of good stage
              col = colg, #good stage in black
              cex = 0.5)
      }
    }
    graphics::legend("bottomleft", legend = c("good", "incipient"), col = c(colg, coli), 
           bty = "n", lwd = 2)
  }
}


plot_tip_traits_df <- function(x){
  
  checkmate::assert_class(x, "list")
  
  if(inherits(x, "cs_result_summaries")){
    df <- data.frame(tip_values = x$tip_traits, 
                     competition = x$competition,
                     selection = x$selection)
    
    return(df)
  }
  
  res_class <- sapply(x, FUN = function(i){inherits(i, "cs_result_summaries")})
  res_class_n <- sum(res_class)
  
  if(sum(!res_class) == length(x)){
    usethis::ui_stop("No objects of class cs_result_summaries in {usethis::ui_field('x')}")
  }
  
  if(res_class_n != length(x)){
    usethis::ui_warn("{res_class_n} of {length(x)} elements of {usethis::ui_field('x')} not objects of class cs_result_summaries. Ignored")
  }
  
  l_df <- lapply(x[res_class], function(i){
    data.frame(tip_values = i$tip_traits, competition = i$competition,
               selection = i$selection)
  })
  return(do.call(rbind, l_df))
}

#' Plot tip trait distribution on competition x selection grid
#'
#' @param x an object of class `cs_result_summaries` or list of `cs_result_summaries` objects
#'
#' @return a plot of tip trait distributions
#' @export
#' @import ggplot2
plot_tip_trait_distribution <- function(x) {
  
  x <- plot_tip_traits_df(x)
  
  x$competition <- factor(x$competition)
  x$selection <- factor(x$selection)
  
  x %>%
    ggplot(aes(x = tip_values)) +
    geom_density(aes(color = competition), alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values =  harrypotter::hp(n = 6, option = "Ravenclaw")) +
    labs(x = "Trait value at tip") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          strip.background.y = element_blank(),
          strip.text.y = element_blank()) +
    facet_grid(rows = vars(competition),
               cols = vars(selection),
               labeller = label_both)
}


plot_trait_tip_means  <- function(x) {
  x <- plot_tip_traits_df(x)
  
  x$competition <- factor(x$competition)
  x$selection <- factor(x$selection)
  
    group_by(competition, selection, iteration) %>%
    summarise(mean = mean(var, na.rm = T),
              sd = stats::sd(var, na.rm = T),
              kurt = e1071::kurtosis(var, na.rm = T)) %>%
    ggplot(aes(y = mean, x = competition, fill = as.factor(competition))) +
    geom_boxplot() +
    geom_smooth(method = "loess", aes(group=1), colour = "dark red", show.legend = F) +
    scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
    labs(x = "Competition", y = "Mean Trait at Tips") +
    #ylab("richness") +
    theme_pubr() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    facet_grid(cols = vars(selection),
               labeller = label_both)
  
  
}

#' Plot result summary variable against time.
#'
#' @param x an object or list of objects of class `cs_result_summaries`.
#' @param variable character string. The variable to be plotted through time.
#' @param ylab character string. Label for y axis. Defaults to variable.
#' @param time character string. If `step`, values are plotted against each time step.
#' If `int` the  values at each integer time point are plotted
#'
#' @return plot of variable against time
#' @export
#' @import ggplot2
plot_var_vs_time <- function(x, variable = c("VAR", "MNND", "VNND"), 
                             ylab = NULL, time = c("step", "int")) {
  
  variable <- match.arg(variable)
  time <- match.arg(time)
  
  if (is.null(ylab)) {
    ylab <- variable
  }
  
  checkmate::assert_class(x, "list")
  
  if(inherits(x, "cs_result_summaries")){
    x <- list(x)
  }
  
  plot_df <- lapply(x, prep_time_plot_df, 
                    variable = variable) %>%
    do.call(rbind, .) 
  
  switch (time,
    step = {
      plot_df$time <- plot_df$time_step
      },
    int = {
      plot_df$time <- plot_df$time_int
      }
  )
  


  plot_df %>%
    dplyr::group_by(.data$selection, 
                    .data$competition,
                    .data$time) %>%
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE)) %>%
    ggplot(aes(x = time, y = value, colour = competition, 
               fill = competition)) +
    geom_line(size = 1) +
    scale_color_manual(values =  harrypotter::hp(n = 6, option = "Ravenclaw"), 
                       name = "Competition") +
    labs(x = "Time", y = ylab) +
      ggpubr::theme_pubr() + 
    facet_grid(rows = vars(selection),
               scales = "free",
               labeller = label_both)
  
}


#' Plot result diversification rates against time.
#'
#' @param x an object or list of objects of class `cs_result_summaries`.
#' @return plot of diversification against time
#' @export
#' @import ggplot2
plot_diversification <- function(x) {
  
  checkmate::assert_class(x, "list")
  
  if(inherits(x, "cs_result_summaries")){
    x <- list(x)
  }
  
  plot_df <- lapply(x, prep_lineages_div) %>%
    mapply(function(x, y) {tibble::add_column(x, sim = y)}, 
           ., seq_along(.), SIMPLIFY = FALSE) %>%
    do.call(rbind, .) %>%
    tidyr::drop_na(.data$spec_ct) 
  
  tree_df <- lapply(x, prep_trees_div) %>%
    mapply(function(x, y) {tibble::add_column(x, sim = y)}, 
           ., seq_along(.), SIMPLIFY = FALSE) %>%
    do.call(rbind, .) 
  
  
  plot_df_spec <- plot_df %>%
    dplyr::mutate(time_bin = as.integer(.data$spec_ct)  + 1L) %>% 
    dplyr::group_by(.data$competition, .data$selection, .data$sim, .data$time_bin) %>%
    dplyr::summarise(rate = dplyr::n()) %>%
    tibble::add_column(process = "speciation") %>%
    dplyr::ungroup()
  
  
  plot_df_ext <- plot_df %>%
  dplyr::filter(.data$status == -2) %>%
    dplyr::mutate(time_bin = as.integer(.data$spec_ct)  + 1L) %>% 
    dplyr::group_by(.data$competition, .data$selection, .data$sim, .data$time_bin) %>%
    dplyr::summarise(rate = dplyr::n()) %>%
    tibble::add_column(process = "extinction")
  
  lin_df_div <- rbind(plot_df_spec, plot_df_ext) %>%
    tidyr::spread(key = .data$process, value = .data$rate, fill = 0) %>%
    dplyr::left_join(tree_df, by = c("sim", "time_bin", 
                                     "competition", "selection")) %>%
    dplyr::mutate(speciation = .data$speciation / .data$branch_len,
           extinction = .data$extinction / .data$branch_len,
           diversification = .data$speciation - .data$extinction,
           turnover = tidyr::replace_na(
             dplyr::na_if(.data$extinction / .data$speciation, Inf), 0)) %>%
    tidyr::gather("process", "rate", c(.data$extinction, .data$speciation, 
                                       .data$diversification)) %>%
    dplyr::group_by(.data$time_bin, .data$competition, .data$selection, .data$process) %>%
    dplyr::summarise(rate = mean(.data$rate)) %>%
    dplyr::ungroup()
  
  lin_df_div %>%
    ggplot(aes(x = as.numeric(time_bin), y = rate, color = as.factor(process))) +
    #geom_smooth(aes(group = as.factor(sim)), color = "light grey", alpha = 0.1, se = F) +
    #geom_point(color = "light grey", size = 0.1, alpha = 0.1) +
    geom_line() +
    # scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50)) +
    scale_color_manual(values = harrypotter::hp(n = 3, option = "DracoMalfoy"), name = "Process") +
    labs(x = "Time", y = "Rate") +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_blank()) +
    facet_grid(rows = vars(competition),
               cols = vars(selection),
               labeller = label_both)
  
  
  
}

# ---- PLOT Data - preprocess ---- ############################################

prep_time_plot_df <- function(x, variable = "VAR") {
  
  df <- tibble::tibble(value = x[[variable]],
                   variable = variable,
                   competition = as.factor(x$competition),
                   selection = as.factor(x$selection))
  
  df$time_step <- seq_along(df$variable) * x$step_size
  df$time_int <- as.integer(df$time_step)
  
  return(df)
}

prep_lineages_div <- function(x) {
  
  df <- as.data.frame(x$lineages)
  df$competition <- x$competition
  df$selection <- x$selection
  
  return(df)
  
}

prep_trees_div <- function(x) {
  
  df <- sumBL(x)
  df$competition <- x$competition
  df$selection <- x$selection
  
  return(df)
  
}
