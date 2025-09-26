library(tidyverse)
library(magrittr)
library(lubridate)
library(grid)  
library(gtable)

get_clusters <- function(snpepi_data, genet_dist, temp_dist, metadata){
    graph <- get_cluster_graph(
        snpepi_data, dist_column = c("dist", "weeks"), "pair_location", 
        dist_threshold = c(genet_dist, temp_dist))
    clusters <- get_cluster_membership_from_graph(graph) %>%
        right_join(metadata %>% select(-any_of("Cluster")), by = "id") 
}

get_only_matching_data <- function(epi_snp_data, metadata, kleborate_data){
    s <- unique(c(epi_snp_data$iso1, epi_snp_data$iso2))
    m <- metadata$id
    matching_ids <- base::intersect(s,m)
    unmatched_ids <- nrow(metadata) - length(matching_ids)
    # if (unmatched_ids > 0){
    print(glue::glue("Excluding {unmatched_ids} genomes without matching SNP data"))
    # }
    mtdt <- metadata %>% filter(id %in% matching_ids)
    kbdt <- kleborate_data %>% filter(`Genome Name` %in% matching_ids)
    snp <- epi_snp_data %>% 
        filter(iso1 %in% matching_ids) %>% 
        filter(iso2 %in% matching_ids) 
    return(list(metadata = mtdt, kleborate_data = kbdt,
                epi_snp_data = snp))
}

get_cluster_details_all_studies <- function(studies_list) {
    studies_list %>% 
        purrr::map_dfr(\(clusters_data) summarise_cluster(clusters_data) %>% 
                           pivot_wider(names_from=name, values_from=value),
                       .id = "Study")
}

plot_sensitivity_SNP_vs_temp_range_2 <- function(
        cluster_and_transmission_sensitivity_df,
        temp_dist_vals=c(1, 2, 4, 8, 52), prop_var='cluster_prop', 
        y_title = "Proportion in clusters", plot_title = NULL, 
        selected_dist_thresh=10, selected_temp_thresh=4,
        range_cols=c("#ffffff", "#8b0000", "#ffc1c1"), 
        facet_col=NULL, facet_rows=1){
    # Use intuitive var names for interactive plot
    y_vars <- paste0(prop_var, " at ", temp_dist_vals, " week(s) threshold")

    df <- cluster_and_transmission_sensitivity_df %>% 
        dplyr::filter(temporal_threshold %in% temp_dist_vals) %>% 
        unique()
    if (is.null(facet_col)){
        df <- df %>% 
            dplyr::select(distance_threshold, temporal_threshold, !!sym(prop_var)) %>% 
            tidyr::pivot_wider(id_cols = c(distance_threshold), 
                               names_from=temporal_threshold, values_from=!!sym(prop_var))
    } else {
        if (!facet_col %in% names(df)) {
            stop(paste("'facet_col' missing: No column named ", facet_col))
        }
        df <- df %>% 
            dplyr::select(distance_threshold, temporal_threshold, !!sym(prop_var), 
                          !!sym(facet_col)) %>% 
            tidyr::pivot_wider(id_cols=c(distance_threshold, !!sym(facet_col)), 
                               names_from=temporal_threshold, values_from=!!sym(prop_var))
    }
    df %<>% dplyr::rename_at(vars(as.character(temp_dist_vals)), list(~y_vars))
    
    # Plot
    p <- df %>% 
        ggplot2::ggplot(aes(x = distance_threshold)) +
        ggplot2::geom_ribbon(aes(ymin = .data[[y_vars[1]]], ymax = .data[[y_vars[5]]]), 
                             fill = range_cols[3]) +
        ggplot2::geom_ribbon(aes(ymin = .data[[y_vars[2]]], ymax = .data[[y_vars[4]]], 
                                 x = distance_threshold), fill = range_cols[2]) +
        ggplot2::geom_vline(aes(xintercept=selected_dist_thresh), lty=2, alpha=.3) +
        ggplot2::geom_line(aes(y = .data[[y_vars[3]]]), colour = range_cols[1]) +
        ggplot2::geom_point(
            data = cluster_and_transmission_sensitivity_df %>% 
                dplyr::filter(distance_threshold == selected_dist_thresh &
                                  temporal_threshold == selected_temp_thresh),
            aes(x=selected_dist_thresh, y=.data[[prop_var]]), colour=range_cols[1], size=2) +
        ggplot2::theme_minimal() + ggplot2::ylim(0, 1) + 
        ggplot2::labs(x = "Genetic distance threshold", y = y_title, title = plot_title) +
        custom_plots_theme 
    if (! is.null(facet_col)){
        p <- p + ggplot2::facet_wrap(facet_col, nrow=facet_rows) + ggplot2::theme_bw()
    }
    return(p)
}


manual_sensitivity_legend <- function(mid_temp_val=4, temp_range_1=c(2,8),
                                      temp_range_2=c(1,52), temp_unit="wks",
                                      legend_title="Estimates at:", l1_type=c("line", "box"),
                                      range_cols=c("#ffffff", "#8b0000", "#ffc1c1")) {
    l1_type = match.arg(l1_type)
    # Grobs
    if (l1_type == "line"){
        L1=rectGrob(height=.08, width=.5, gp=gpar(fill=range_cols[1], col=range_cols[2], lwd=0.5))
    } else if (l1_type == "box") {
        L1=rectGrob(height=.5, width=.5, gp=gpar(fill=range_cols[1], col=range_cols[2]))
    }
    L2=rectGrob(height=.5, width=.5, gp=gpar(fill=range_cols[2], col=NA))
    L3=rectGrob(height=.5, width=.5, gp=gpar(fill=range_cols[3], col=NA))
    T1=textGrob(paste(mid_temp_val, temp_unit),
                x=.2, just="left", gp=gpar(fontsize=12, lineheight=0.9))
    T2=textGrob(paste(paste0(temp_range_1, collapse=" - "), temp_unit),
                x=.2, just="left", gp=gpar(fontsize=12, lineheight=0.9))
    T3=textGrob(paste(paste0(temp_range_2, collapse=" - "), temp_unit),
                x=.2, just="left", gp=gpar(fontsize=12, lineheight=0.9))
    # Construct gtable - 2 columns X 4 rows
    leg=gtable(width=unit(c(1,3), "cm"), height=unit(c(1.2,1.2,1.2,1.2), "cm"))
    leg=gtable_add_grob(leg, rectGrob(gp=gpar(fill="white", col=NA)), t=2,l=1,b=4,r=2)
    # Place the six grobs into the table
    leg=gtable_add_grob(leg, list(L1,L2,L3,T1,T2,T3),
                        t=c(2,3,4,2,3,4), l=c(1,1,1,2,2,2))
    # Title
    if (!is.null(legend_title)) {
        leg=gtable_add_grob(leg, textGrob(legend_title, gp=gpar(fontsize=14)),
                            t=1, l=1, r=2)
    }
    return(leg)
}

manual_sensitivity_legend2 <- function(mid_temp_val=4, temp_range_1=c(2,8), temp_range_2=c(1,52), 
                                       selected_dist_thresh=10, selected_temp_thresh=4,
                                       dist_unit="SNVs", temp_unit="wks",
                                       legend_title="Estimates at: ", l1_type=c("line", "box"),
                                       range_cols=c("#ffffff", "#8b0000", "#ffc1c1"), fontsize=10) {
    library(grid)
    library(gtable)
    l1_type = match.arg(l1_type)
    # Grobs
    if (l1_type == "line"){
        L1=rectGrob(height=.06, width=.5, gp=gpar(fill=range_cols[1], col=range_cols[2], lwd=0.5))
    } else if (l1_type == "box") {
        L1=rectGrob(height=.5, width=.5, gp=gpar(fill=range_cols[1], col=range_cols[2]))
    }
    L2=rectGrob(height=.5, width=.5, gp=gpar(fill=range_cols[2], col=NA))
    L3=rectGrob(height=.5, width=.5, gp=gpar(fill=range_cols[3], col=NA))
    L4=circleGrob(x=.5, y=.5, r=.25, gp=gpar(fill=range_cols[1], col=range_cols[2], lwd=0.4))
    
    T1=textGrob(paste(mid_temp_val, temp_unit),
                x=.2, just="left", gp=gpar(fontsize=fontsize, lineheight=0.8))
    T2=textGrob(paste(paste0(temp_range_1, collapse=" - "), temp_unit),
                x=.2, just="left", gp=gpar(fontsize=fontsize, lineheight=0.8))
    T3=textGrob(paste(paste0(temp_range_2, collapse=" - "), temp_unit),
                x=.2, just="left", gp=gpar(fontsize=fontsize, lineheight=0.8))
    T4=textGrob(glue::glue("{selected_temp_thresh} {temp_unit} &\n {selected_dist_thresh} {dist_unit}"),
                x=.2, just="left", gp=gpar(fontsize=fontsize, lineheight=0.8))
    # Construct gtable - 2 columns X 5 rows
    leg=gtable(width=unit(c(1,2), "cm"), height=unit(c(1.2,1.2,1.2,1.2,1.2), "cm"))
    leg=gtable_add_grob(leg, rectGrob(gp=gpar(fill="white", col=NA)), t=2,l=1,b=5,r=2)
    # Place the six grobs into the table
    leg=gtable_add_grob(leg, list(L1,L2,L3,L4, T1,T2,T3,T4),
                        t=c(2,3,4,5,2,3,4,5), l=c(1,1,1,1,2,2,2,2))
    # Title
    if (!is.null(legend_title)) {
        leg=gtable_add_grob(leg, textGrob(legend_title, gp=gpar(fontsize=fontsize)),
                            t=1, l=1, r=2)
    }
    return(leg)
}
    
plot_dist_method_sens_comp <- function(
        d, snp_val, temp_dist_vals, comparison_var="comparison_group",
        prop_var='cluster_prop', y_title="Proportion in clusters", 
        range_cols=c("white", "#8b0000", "#ffc1c1"),
        plot_title=NULL){
    stopifnot(length(snp_val) == 1, length(temp_dist_vals) == 5)
    d %<>% rename("Group" := !!rlang::sym(comparison_var))
    # arrange sites by prop_var at midpoint temp_dist
    ordered_d <- d %>% 
        filter(temporal_threshold == temp_dist_vals[3]) %>% 
        arrange(study, desc(.data[[prop_var]]))
    # Use intuitive var names for interactive plot
    y_vars <- paste0(prop_var, " at ", temp_dist_vals, " week(s) threshold") 
    # plot data
    d <- d %>% 
        select(distance_threshold, temporal_threshold, dist_method, Group, all_of(c(prop_var))) %>% 
        mutate(Group=factor(Group, levels=ordered_d$Group)) %>% 
        pivot_wider(id_cols=c(distance_threshold, dist_method, Group), 
                    names_from=temporal_threshold, 
                    values_from=all_of(c(prop_var))) %>% 
        rename_at(vars(as.character(temp_dist_vals)), list(~y_vars))
    
    # plot
    p <- d %>% 
        ggplot(aes(x=Group)) +
        geom_linerange(aes(ymin=.data[[y_vars[1]]], ymax=.data[[y_vars[5]]]), 
                       col=range_cols[3], lwd=4) +
        geom_linerange(aes(ymin=.data[[y_vars[2]]], ymax=.data[[y_vars[4]]]), 
                       col=range_cols[2], lwd=2.5) +
        geom_point(aes(y=.data[[y_vars[3]]]), shape=22, color=range_cols[2], fill=range_cols[1], size=2) +
        geom_text(aes(label=glue::glue("{.data[[y_vars[3]]]} (2-8 wks: {.data[[y_vars[2]]]}-{.data[[y_vars[4]]]})
                                       (1-52 wks: {.data[[y_vars[1]]]}-{.data[[y_vars[5]]]})"), 
                      y=.data[[y_vars[3]]]), 
                  size=4, hjust=0.5, vjust=1.4) 
    # Plot SKA estimates with diff colours
    d_ska <- d %>% filter(dist_method == "SKA")
    if (nrow(d_ska) > 0) {
        p <- p + 
            geom_linerange(data=d_ska, aes(ymin=.data[[y_vars[1]]], ymax=.data[[y_vars[5]]]),
                           col=range_cols[3], lwd=4) +
            geom_linerange(data=d_ska, aes(ymin=.data[[y_vars[2]]], ymax=.data[[y_vars[4]]]), 
                           col=range_cols[2], lwd=2.5) +
            geom_point(data=d_ska, aes(y=.data[[y_vars[3]]]), 
                       shape=22, color=range_cols[2], fill=range_cols[1], size=2)
    }
    # Finish plot
    p <- p + theme_minimal() + ylim(0, 1) + 
        labs(x=NULL, y=y_title) +
        custom_plots_theme +
        theme(plot.title=element_text(face='bold', size=14),
              axis.text=element_text(size=12),
              axis.title=element_text(size=14)) +
        coord_flip() +
        ggtitle(plot_title)
    
    return(p)
}

summarise_cluster_full <- function(clusters_data) {
    clust_info <- get_cluster_info(clusters_data)
    tot_isolates <- n_distinct(clusters_data$id)
    n_clusters <- n_distinct(clusters_data$Cluster, na.rm = T)
    n_iso_in_clust <- clusters_data %>% filter(!is.na(Cluster)) %>% nrow()
    n_iso_transm <- n_iso_in_clust - n_clusters
    clust_size <- clust_info %>% 
        reframe(x = paste0(
            "Median=",median(`N isolates`,na.rm = T), "; ",
            "Range=", min(`N isolates`,na.rm = T), "-", 
            max(`N isolates`,na.rm = T))) %>% 
        pull(x)
    clust_duration <- clust_info %>% 
        reframe(x = paste0(
            "Median=", median(`Duration (days)`,na.rm = T), "; ",
            "Range=", min(`Duration (days)`,na.rm = T), "-", 
            max(`Duration (days)`,na.rm = T))) %>% 
        pull(x)
    n_clust_st <- clusters_data %>% dplyr::filter(!is.na(Cluster)) %>% 
        pull(ST) %>% n_distinct(na.rm = T)
    
    summary <- tibble::tribble(
        ~name, ~value,
        'Total isolates', as.character(tot_isolates),
        'N clusters', as.character(n_clusters),
        'N isolates in clusters', as.character(n_iso_in_clust),
        'N transmitted isolates', as.character(n_iso_transm),
        'Cluster size', as.character(clust_size),
        'Duration (days)', as.character(clust_duration),
        'Cluster STs', as.character(n_clust_st)
    )
    return(summary)
}

summarise_cluster_nodates <- function(clusters_data) {
    tot_isolates <- n_distinct(clusters_data$id)
    n_clusters <- n_distinct(clusters_data$Cluster, na.rm = T)
    n_iso_in_clust <- clusters_data %>% filter(!is.na(Cluster)) %>% nrow()
    n_iso_transm <- n_iso_in_clust - n_clusters
    clust_size <- clusters_data %>% filter(!is.na(Cluster)) %>% group_by(Cluster) %>% 
        reframe("N isolates" = n()) %>% 
        reframe(x = paste0(
            "Median=",median(`N isolates`,na.rm = T), "; ",
            "Range=", min(`N isolates`,na.rm = T), "-", 
            max(`N isolates`,na.rm = T))) %>% 
        pull(x)
    summary <- tibble::tribble(
        ~name, ~value,
        'Total isolates', as.character(tot_isolates),
        'Clusters', as.character(n_clusters),
        'N isolates in clusters', as.character(n_iso_in_clust),
        'N transmitted isolates', as.character(n_iso_transm),
        'Cluster size', as.character(clust_size),
        '# Cluster STs', clusters_data %>% dplyr::filter(!is.na(Cluster)) %>% 
            pull(ST) %>% n_distinct(na.rm = T) %>% as.character()
    )
    return(summary)
}

get_cluster_info_nodates <- function(clusters_data){
    clusters_data %>% 
        dplyr::filter(!is.na(Cluster)) %>% 
        dplyr::group_by(Cluster) %>% 
        dplyr::reframe(Site = paste0(sort(unique(Site)), collapse = '; '),
                       "N isolates" = n(),
                       ST = paste0(sort(unique(ST)), collapse = '; ')) %>% 
        dplyr::arrange(desc(`N isolates`))
}

summarise_ko_by_case_type <- function(deduplicated_clusters_data, column_to_summarise){
    deduplicated_clusters_data %<>% 
        rename("col_to_summarise" := !!sym(column_to_summarise))
    df <- deduplicated_clusters_data %>% 
        count(col_to_summarise) %>% 
        mutate(pct = round(n / sum(n) * 100,1)) %>% 
        mutate(All = paste0(n, " (", pct, "%)")) %>% 
        arrange(desc(pct)) %>% select(col_to_summarise, All) %>% 
        left_join(
            deduplicated_clusters_data %>% 
                mutate(case_type = if_else(is.na(Cluster), "Singletons", "Clusters")) %>% 
                group_by(case_type, col_to_summarise) %>% 
                summarise(N = n()) %>% 
                mutate(total_isolates = sum(N, na.rm = TRUE),  
                       pct = round((N / total_isolates) * 100, 1)) %>% 
                select(-total_isolates) %>% ungroup() %>%  
                mutate(v = paste0(N, " (", pct, "%)")) %>% 
                pivot_wider(id_cols = col_to_summarise, names_from=case_type, values_from=v),
            by = "col_to_summarise") %>% 
        left_join(
            deduplicated_clusters_data %>% 
                group_by(col_to_summarise) %>% 
                reframe(
                    n_countries = n_distinct(Country, na.rm=T),
                    n_countries_clusters = n_distinct(Country[!is.na(Cluster)], na.rm=T),
                    n_countries_singletons = n_distinct(Country[is.na(Cluster)], na.rm=T),
                    n_sites = n_distinct(Site, na.rm=T),
                    n_sites_clusters = n_distinct(Site[!is.na(Cluster)], na.rm=T),
                    n_sites_singletons = n_distinct(Site[is.na(Cluster)], na.rm=T),
                    n_clones = n_distinct(ST, na.rm=T),
                    n_clones_clusters = n_distinct(ST[!is.na(Cluster)], na.rm=T),
                    n_clones_singletons = n_distinct(ST[is.na(Cluster)], na.rm=T)
                ),
            by = "col_to_summarise"
        )
    df %>% rename(!!sym(column_to_summarise) := "col_to_summarise")
}

plot_facility_data_heatmap <- function(facilities_df, meta_df){
    facs_df <- facilities_df %>% 
        select(Site, Location, Facility_type, Facility_size, 
               Neonatal_beds, Piped_water, Onsite_surgery) %>%
        mutate(Site = factor(Site, levels = meta_df$S))
    # colors
    facs_cols_loc <- c("Rural"="#e5cb9a", "Urban"="#000000")
    facs_cols_surg <- c("No"="#e5cb9a", "Yes"="#000000")
    # facs_cols_type <- c("District"="#cfe2f3", "Regional"="#2986cc", "Tertiary"="#1b3b5a")
    facs_cols_type <- c("District"="#e5cb9a", "Regional"="#7a2828", "Tertiary"="#000000")
    facs_cols_size <- c("Small" = "#e5cb9a", "Medium" = "#dda960", 
                        "Large" = "#7a2828", "Very Large" = "#000000")
    facs_cols_water <- c("Sometimes" = "#e5cb9a", "Most times" = "#7a2828", "Always" = "#000000")
    # theme
    facs_plots_theme <- function(p) {
        p + coord_cartesian(ylim=c(-1,nrow(meta_df)+2)) + 
            scale_x_discrete(expand=c(0,0)) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle=45, hjust=1, size=14, lineheight=0.75), 
                  legend.text = element_text(size=14), legend.title = element_text(size=14),
                  axis.title=element_blank(), axis.text.y=element_blank(),
                  legend.key.height = unit(2, "mm"),
                  plot.margin = unit(c(0.2, 0.2, 0.2, 0.05), "cm"))
    }
    # Plots ----------------------------
    # Location
    p_loc <- facs_df %>% select(Site, Location) %>% 
        pivot_longer(-Site) %>% 
        ggplot(aes(y=Site, x=name, fill=value)) + geom_tile() +
        scale_fill_manual(values=facs_cols_loc, name="Location", na.value="white")
    p_loc <- facs_plots_theme(p_loc)
    # Type
    p_type <- facs_df %>% select(Site, Type = Facility_type) %>% 
        pivot_longer(-Site) %>% 
        ggplot(aes(y=Site, x=name, fill=value)) + geom_tile() +
        scale_fill_manual(values=facs_cols_type, name="Facility type", na.value="white")
    p_type <- facs_plots_theme(p_type)
    # Beds
    p_beds <- facs_df %>% select(Site, `No. beds` = Neonatal_beds) %>% 
        pivot_longer(-Site) %>% 
        ggplot(aes(y=Site, x=name, fill=value)) + geom_tile() +
        scale_fill_viridis_c(
            option='B', name="No. neonatal beds", oob = scales::squish,
            na.value="white")
    p_beds <- facs_plots_theme(p_beds)
    # Size
    p_size <- facs_df %>% select(Site, Size = Facility_size) %>% 
        pivot_longer(-Site) %>% 
        ggplot(aes(y=Site, x=name, fill=value)) + geom_tile() +
        scale_fill_manual(values=facs_cols_size, name="Facility size", na.value="white")
    p_size <- facs_plots_theme(p_size) 
    # Piped water 
    p_water <- facs_df %>% select(Site, `Piped water` = Piped_water) %>% 
        pivot_longer(-Site) %>% 
        ggplot(aes(y=Site, x=name, fill=value)) + geom_tile() +
        scale_fill_manual(values=facs_cols_water, name="Piped water available",
                          na.value="white")
    p_water <- facs_plots_theme(p_water)
    # Onsite surgery 
    p_surg <- facs_df %>% select(Site, `Onsite neonatal\nsurgical facilities` = 
                                     Onsite_surgery) %>% 
        pivot_longer(-Site) %>% 
        ggplot(aes(y=Site, x=name, fill=value)) + geom_tile() +
        scale_fill_manual(values=facs_cols_surg, na.value="white",
                          name="Onsite neonatal surgical facilities")
    p_surg <- facs_plots_theme(p_surg)
    
    return(list(p1=p_loc, p2=p_type, p3=p_surg, p4=p_size, p5=p_beds, p6=p_water))
}

plot_clust_prop_by_facility_data <- function(clust_prop_and_facilities_df, vars_to_plot = c("Location")){
    my_plots_list <- list()
    for (fac in vars_to_plot) {
        lab = gsub("_"," ", fac)
        p <- clust_prop_and_facilities_df %>% 
            ggplot(aes(y = cluster_prop, x = .data[[fac]])) +
            geom_boxplot() + geom_jitter(width=0.1) +
            scale_y_continuous(limits = c(0,1)) +
            labs(y = "Cluster proportion", x = NULL, title=lab) +
            theme_bw() +
            theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
                  plot.title=element_text(size=14))
        my_plots_list[[fac]] <- p
    }
    return(my_plots_list)
}


parse_logistf_model <- function(object){
    # See getAnywhere("summary.logistf")
    if (!is.null(object$modcontrol$terms.fit)) {
        var.red <- object$var[object$modcontrol$terms.fit, object$modcontrol$terms.fit]
        coefs <- coef(object)[object$modcontrol$terms.fit]
        chi2 <- vector(length = length(object$terms))
        chi2[object$modcontrol$terms.fit] <- 
            qchisq(1 - object$prob[object$modcontrol$terms.fit], 1)
        chi2[-object$modcontrol$terms.fit] <- 0
    }
    else {
        var.red <- object$var
        coefs <- coef(object)
        chi2 <- qchisq(1 - object$prob, 1)
    }
    out <- cbind(object$coefficients, diag(object$var)^0.5, object$ci.lower, 
                 object$ci.upper, chi2, object$prob, 
                 ifelse(object$method.ci == "Wald", 1, ifelse(object$method.ci == "-", 3, 2)))
    dimnames(out) <- list(names(object$coefficients), 
                          c("Estimate", "StdErr", 
                            paste0(c("Lower_", "Upper_"), 1 - object$alpha),
                            "Chisq", "P", "Method"))
    d <- out[-1,] %>% 
        as_tibble(rownames = "Predictor") %>% 
        mutate(OR = round(exp(Estimate), 2),
               lower = round(exp(Estimate - (1.96 * StdErr)), 2),
               upper = round(exp(Estimate + (1.96 * StdErr)), 2)
        ) %>% 
        mutate(P = case_when(P > 0.05 ~ round(P, 3), 
                             P > 0.0001 ~ round(P, 5), 
                             TRUE ~ P)) %>% 
        select(Predictor, Estimate, StdErr, OR, lower, upper, P)
    return(d)
}

parse_metareg_model <- function(object){
    out <- cbind(object$beta, object$se, object$pval)
    colnames(out) <- c("Estimate", "StdErr", "P")
    d <- out[-1,] %>% 
        as_tibble(rownames = "Predictor") %>% 
        mutate(OR = round(exp(Estimate), 2),
               lower = round(exp(Estimate - (1.96 * StdErr)), 2),
               upper = round(exp(Estimate + (1.96 * StdErr)), 2)
        ) %>% 
        mutate(P = case_when(P > 0.05 ~ round(P, 3), 
                             P > 0.0001 ~ round(P, 5), 
                             TRUE ~ P)) %>% 
        select(Predictor, Estimate, StdErr, OR, lower, upper, P)
    return(d)
}

get_lm_eqn <- function(df, format=F){
    m <- lm(y ~ x -1, df);
    slope <- unname(coef(m)[1])
    r_squared <- summary(m)$r.squared
    if(format){
        # Use in plot with geom_richtext or equivalent fnx
        return(sprintf("Slope = %.2f, R<sup>2</sup> = %.3f", slope, r_squared))
    } else {
        return(sprintf("Slope = %.2f, R2 = %.3f", slope, r_squared))
    }
}
