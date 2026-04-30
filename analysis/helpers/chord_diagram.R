#!/usr/bin/env Rscript
# chord_diagram.R — reusable circular chord/PPI network diagram
#
# Extracted from figure1.R panel F. Call make_chord_ppi() to produce a ggplot
# circular layout identical in style to the pre-onset DEP network in Figure 1.
#
# Dependencies: igraph, data.table, ggplot2 (all assumed loaded by caller)

# ── Bezier helper (cubic) ─────────────────────────────────────────────────────
.bzier <- function(p0, p1, p2, p3, n = 80) {
    t <- seq(0, 1, length.out = n)
    data.frame(
        x = (1-t)^3*p0[1] + 3*(1-t)^2*t*p1[1] + 3*(1-t)*t^2*p2[1] + t^3*p3[1],
        y = (1-t)^3*p0[2] + 3*(1-t)^2*t*p1[2] + 3*(1-t)*t^2*p2[2] + t^3*p3[2]
    )
}

# ── make_chord_ppi() ──────────────────────────────────────────────────────────
#
# Arguments:
#   g              igraph — undirected, edge attribute `combined_score` (0–1)
#   comm_vec       named integer: node name → community ID (1-based)
#                  community 0 = orphan (not shown in sectors)
#   comm_cols      named character: community_id → sector/chord colour
#   comm_cols_light named character: community_id → light fill for node circles
#   comm_labels    named character: community_id → text label
#   node_direction named character: node name → "Up" / "Down" / "NS"
#                  determines node stroke colour (Up = warm, Down = cool, NS = grey)
#   max_nodes      integer — prune to top N by degree before plotting
#   R_RING         ring radius (default 3.0)
#   GAP            gap between sectors as fraction of 2π (default 0.012)
#   title_str      plot title
#   subtitle_str   plot subtitle

make_chord_ppi <- function(g,
                            comm_vec,
                            comm_cols,
                            comm_cols_light,
                            comm_labels,
                            node_direction = NULL,
                            max_nodes      = 60L,
                            R_RING         = 3.0,
                            GAP            = 0.012,
                            title_str      = "",
                            subtitle_str   = "") {

    # Prune to top N by degree
    if (vcount(g) > max_nodes) {
        deg_all <- degree(g)
        keep    <- names(sort(deg_all, decreasing = TRUE))[seq_len(max_nodes)]
        g       <- induced_subgraph(g, keep)
    }
    deg <- degree(g)

    # Resolve communities for nodes actually in the graph
    cid <- comm_vec[V(g)$name]
    cid[is.na(cid)] <- 0L

    # Direction for node stroke
    if (is.null(node_direction)) {
        dir_map <- setNames(rep("NS", vcount(g)), V(g)$name)
    } else {
        dir_map <- node_direction[V(g)$name]
        dir_map[is.na(dir_map)] <- "NS"
    }

    n_c <- max(cid[cid > 0])
    if (n_c < 1L) stop("No community IDs > 0 found")

    # ── Radial positions ──────────────────────────────────────────────────────
    node_dt <- data.table(
        name      = V(g)$name,
        comm      = cid,
        degree    = deg,
        direction = dir_map
    )
    setorder(node_dt, comm, -degree)

    N_top  <- sum(node_dt$comm > 0)
    avail  <- 1.0 - GAP * n_c
    pos    <- vector("list", nrow(node_dt)); names(pos) <- node_dt$name
    cang   <- vector("list", n_c)
    cur    <- 0.0

    for (ci in seq_len(n_c)) {
        ii <- which(node_dt$comm == ci); n_ci <- length(ii)
        if (n_ci == 0L) next
        frac <- avail * n_ci / N_top
        sa <- cur * 2*pi - pi/2; ea <- (cur + frac) * 2*pi - pi/2
        for (j in seq_len(n_ci)) {
            a <- (cur + frac * (j - 0.5) / n_ci) * 2*pi - pi/2
            pos[[node_dt$name[ii[j]]]] <- c(R_RING * cos(a), R_RING * sin(a))
        }
        cang[[ci]] <- list(s = sa, e = ea, m = (sa + ea) / 2)
        cur <- cur + frac + GAP
    }

    node_dt[, x := vapply(pos[name], function(p) if (is.null(p)) NA_real_ else p[1], numeric(1))]
    node_dt[, y := vapply(pos[name], function(p) if (is.null(p)) NA_real_ else p[2], numeric(1))]
    node_dt <- node_dt[!is.na(x)]
    node_dt[, comm_col       := comm_cols[as.character(comm)]]
    node_dt[, comm_col_light := comm_cols_light[as.character(comm)]]

    # ── Sector polygons ───────────────────────────────────────────────────────
    sector_list <- lapply(seq_len(n_c), function(ci) {
        if (is.null(cang[[ci]])) return(NULL)
        th <- seq(cang[[ci]]$s, cang[[ci]]$e, length.out = 55)
        lo <- R_RING - 0.55; hi <- R_RING + 0.55
        data.table(x = c(hi*cos(th), lo*cos(rev(th))),
                   y = c(hi*sin(th), lo*sin(rev(th))), ci = ci)
    })

    # ── Edges ─────────────────────────────────────────────────────────────────
    el  <- as_edgelist(g)
    edt <- data.table(u = el[,1], v = el[,2], s = E(g)$combined_score)
    edt[, cu := cid[u]]; edt[, cv := cid[v]]
    edt[, within := (cu == cv & cu > 0)]
    slo <- min(edt$s); shi <- max(edt$s); sR <- max(shi - slo, 1e-3)

    intra_list <- lapply(seq_len(n_c), function(ci) {
        sub <- edt[within == TRUE & cu == ci]; if (nrow(sub) == 0) return(NULL)
        rbindlist(lapply(seq_len(nrow(sub)), function(i) {
            pu <- pos[[sub$u[i]]]; pv <- pos[[sub$v[i]]]
            if (is.null(pu) || is.null(pv)) return(NULL)
            data.table(x = pu[1], y = pu[2], xend = pv[1], yend = pv[2],
                       lw = 0.5 + 1.8 * (sub$s[i] - slo) / sR, ci = ci)
        }))
    })

    inter_pairs <- edt[within == FALSE & cu > 0 & cv > 0,
                       .(tot = sum(s)), by = .(c1 = pmin(cu, cv), c2 = pmax(cu, cv))]
    mb <- max(inter_pairs$tot, 1); Ra <- R_RING - 0.35

    chord_list <- lapply(seq_len(nrow(inter_pairs)), function(k) {
        c1 <- inter_pairs$c1[k]; c2 <- inter_pairs$c2[k]
        if (is.null(cang[[c1]]) || is.null(cang[[c2]])) return(NULL)
        m1 <- cang[[c1]]$m; m2 <- cang[[c2]]$m
        x1 <- Ra*cos(m1); y1 <- Ra*sin(m1); x2 <- Ra*cos(m2); y2 <- Ra*sin(m2)
        cx <- (x1+x2)*0.08; cy <- (y1+y2)*0.08; dx <- x2-x1; dy <- y2-y1
        bz <- .bzier(c(x1,y1), c(cx-dy*0.04, cy+dx*0.04),
                     c(cx+dy*0.04, cy-dx*0.04), c(x2,y2))
        bz$nb <- inter_pairs$tot[k] / mb; bz$c1 <- c1; bz$pid <- k; bz
    })

    # ── Community labels ──────────────────────────────────────────────────────
    R_LBL <- R_RING + 1.20
    clab_dt <- rbindlist(lapply(seq_len(n_c), function(ci) {
        if (is.null(cang[[ci]])) return(NULL)
        m  <- cang[[ci]]$m; lx <- R_LBL*cos(m); ly <- R_LBL*sin(m)
        data.table(x = lx, y = ly, label = comm_labels[as.character(ci)],
                   hjust = if (cos(m) > 0.15) 0 else if (cos(m) < -0.15) 1 else 0.5,
                   col = comm_cols[as.character(ci)])
    }))

    # ── Assemble ggplot ───────────────────────────────────────────────────────
    th_f   <- seq(0, 2*pi, length.out = 300)
    ring_d <- data.frame(x = R_RING*cos(th_f), y = R_RING*sin(th_f))

    p <- ggplot() +
        coord_equal(clip = "off") +
        theme_void() +
        theme(
            plot.title    = element_text(size = 9),
            plot.subtitle = element_text(size = 7, colour = "grey40"),
            legend.position = "none",
            plot.margin   = margin(40, 90, 40, 60, "pt")
        )

    for (ci in seq_len(n_c)) {
        sec <- sector_list[[ci]]; if (is.null(sec)) next
        p <- p + geom_polygon(data = sec, aes(x=x, y=y),
                              fill = comm_cols[as.character(ci)], colour = NA, alpha = 0.30)
    }
    p <- p + geom_path(data = ring_d, aes(x=x, y=y),
                       colour = "#d8d8d8", linewidth = 0.35, linetype = "dashed")

    for (k in seq_along(chord_list)) {
        bz <- chord_list[[k]]; if (is.null(bz)) next
        p <- p + geom_path(data = bz, aes(x=x, y=y),
                           colour = comm_cols[as.character(bz$c1[1])],
                           linewidth = 1.0 + 4.0*bz$nb[1],
                           alpha = 0.20 + 0.22*bz$nb[1])
    }
    for (ci in seq_len(n_c)) {
        seg <- intra_list[[ci]]; if (is.null(seg)) next
        p <- p + geom_segment(data = seg, aes(x=x, y=y, xend=xend, yend=yend),
                              colour = comm_cols[as.character(ci)],
                              linewidth = 0.85, alpha = 0.55)
    }
    p <- p +
        geom_point(data = node_dt, aes(x=x, y=y, fill=comm_col_light),
                   shape = 21, size = 6.5, stroke = 0.65,
                   colour = "grey20", alpha = 0.95) +
        scale_fill_identity() +
        geom_text(data = node_dt, aes(x=x, y=y, label=name),
                  size = 1.55, colour = "grey10", fontface = "bold")

    for (i in seq_len(nrow(clab_dt))) {
        p <- p + annotate("text",
                          x = clab_dt$x[i], y = clab_dt$y[i],
                          label  = clab_dt$label[i],
                          colour = clab_dt$col[i],
                          hjust  = clab_dt$hjust[i],
                          size = 2.4, fontface = "bold", lineheight = 0.85)
    }
    p + labs(title = title_str, subtitle = subtitle_str)
}
