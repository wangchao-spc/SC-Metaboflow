"""
R Bridge Module — Python interface to R packages via Rscript subprocess.

Calls Rscript.exe directly (no rpy2 needed — avoids "R not built as a library"
issue on Windows). Provides wrappers for Seurat, Monocle, EigenMS.
"""

import os
import sys
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# ============================================================
# Find R installation
# ============================================================

_R_HOME = None
_RSCRIPT = None

_POSSIBLE_R_PATHS = [
    r"C:\Program Files\R\R-4.6.1",
    r"C:\Program Files\R\R-4.4.2",
    r"C:\Program Files\R\R-4.4.1",
    r"C:\Program Files\R\R-4.4.0",
    r"C:\Program Files\R\R-4.3.3",
    "/usr/lib/R",
    "/usr/local/lib/R",
]

# Scan Program Files\R for installed versions
_r_base = r"C:\Program Files\R"
if os.path.isdir(_r_base):
    for _entry in sorted(os.listdir(_r_base), reverse=True):
        _full = os.path.join(_r_base, _entry)
        if os.path.isdir(_full) and _entry.startswith("R-"):
            _POSSIBLE_R_PATHS.insert(0, _full)

for _rp in _POSSIBLE_R_PATHS:
    _rscript = os.path.join(_rp, "bin", "x64", "Rscript.exe")
    if os.path.exists(_rscript):
        _R_HOME = _rp
        _RSCRIPT = _rscript
        break
    _rscript = os.path.join(_rp, "bin", "Rscript.exe")
    if os.path.exists(_rscript):
        _R_HOME = _rp
        _RSCRIPT = _rscript
        break

_R_AVAILABLE = _RSCRIPT is not None


def _get_rscript() -> Optional[str]:
    """Get path to Rscript executable."""
    return _RSCRIPT


# ============================================================
# R Bridge Class
# ============================================================


class RBridge:
    """Bridge to R environment using Rscript subprocess."""

    def __init__(self):
        self._r_log = []
        self._seurat_ok = False
        self._monocle_ok = False
        self._eigenms_ok = False

    @property
    def is_available(self) -> bool:
        return _R_AVAILABLE

    @property
    def is_ready(self) -> bool:
        return any([self._seurat_ok, self._monocle_ok, self._eigenms_ok])

    @property
    def status_message(self) -> str:
        if not _R_AVAILABLE:
            return "R not found. Please install R from https://cran.r-project.org/"
        if not self.is_ready:
            return (
                "R found but packages not loaded. "
                "Click 'Install R Packages' to install Seurat/monocle/EigenMS."
            )
        pkgs = []
        if self._seurat_ok:
            pkgs.append("Seurat")
        if self._monocle_ok:
            pkgs.append("Monocle")
        if self._eigenms_ok:
            pkgs.append("EigenMS")
        else:
            pkgs.append("EigenMS (built-in SVD)")
        return f"R Ready — {', '.join(pkgs)} available"

    @property
    def package_status(self) -> dict:
        return {
            "Seurat": self._seurat_ok,
            "Monocle": self._monocle_ok,
            "EigenMS": self._eigenms_ok,
        }

    def get_log(self) -> str:
        return "\n".join(self._r_log)

    def _r_print(self, msg: str):
        self._r_log.append(msg)
        logger.info(msg)

    def _run_r_script(self, script: str, timeout: int = 1800) -> tuple[int, str, str]:
        """Run an R script via Rscript and return (returncode, stdout, stderr)."""
        if not _RSCRIPT:
            return -1, "", "Rscript not found"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".R", delete=False, encoding="utf-8"
        ) as f:
            f.write(script)
            tmp_path = f.name
        try:
            proc = subprocess.run(
                [_RSCRIPT, "--no-save", "--no-restore", tmp_path],
                capture_output=True, text=True, timeout=timeout,
                env={**os.environ, "R_HOME": _R_HOME},
            )
            return proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired:
            return -1, "", "TIMEOUT after {} seconds".format(timeout)
        finally:
            try:
                os.unlink(tmp_path)
            except Exception:
                pass

    # ============================================================
    # R Package Installation
    # ============================================================

    def ensure_r_packages(self) -> bool:
        """Install and verify R packages."""
        if not _R_AVAILABLE:
            self._r_print("ERROR: R not found")
            return False

        self._r_print("=" * 50)
        self._r_print("Installing R packages via Rscript...")
        self._r_print(f"Rscript: {_RSCRIPT}")
        self._r_print("=" * 50)

        # Step 1: Install remotes + BiocManager
        self._r_print("\n[Step 1] Installing CRAN tools...")
        rc, out, err = self._run_r_script('''
            options(repos = c(CRAN = "https://cran.r-project.org"))
            options(install.packages.compile.from.source = "never")

            cat("Installing remotes...\\n")
            if (!requireNamespace("remotes", quietly = TRUE)) {
                install.packages("remotes")
            }
            library(remotes)
            cat("remotes: OK\\n")

            cat("Installing BiocManager...\\n")
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
            }
            library(BiocManager)
            cat("BiocManager: OK\\n")

            cat("Step 1 complete\\n")
        ''', timeout=600)
        self._r_print(out[-500:] if len(out) > 500 else out)
        if err:
            self._r_print(f"STDERR: {err[-300:]}")

        # Step 2: Install CRAN packages
        self._r_print("\n[Step 2] Installing CRAN packages...")
        rc2, out2, err2 = self._run_r_script('''
            options(repos = c(CRAN = "https://cran.r-project.org"))
            options(install.packages.compile.from.source = "never")
            for (pkg in c("ggplot2", "dplyr", "Matrix")) {
                cat("Checking", pkg, "...\\n")
                tryCatch({
                    library(pkg, character.only = TRUE)
                    cat(pkg, ": already installed\\n")
                }, error = function(e) {
                    cat("Installing", pkg, "...\\n")
                    install.packages(pkg)
                    library(pkg, character.only = TRUE)
                    cat(pkg, ": OK\\n")
                })
            }
            cat("Step 2 complete\\n")
        ''', timeout=600)
        self._r_print(out2[-500:] if len(out2) > 500 else out2)
        if err2:
            self._r_print(f"STDERR: {err2[-300:]}")

        # Step 3: Install Bioconductor packages
        self._r_print("\n[Step 3] Installing Bioconductor packages (Seurat, monocle)...")
        rc3, out3, err3 = self._run_r_script('''
            options(repos = c(CRAN = "https://cran.r-project.org"))
            options(install.packages.compile.from.source = "never")
            library(BiocManager)

            for (pkg in c("monocle", "Seurat")) {
                cat("Checking", pkg, "...\\n")
                tryCatch({
                    library(pkg, character.only = TRUE)
                    cat(pkg, ": already installed\\n")
                }, error = function(e) {
                    cat("Installing", pkg, "(this may take 5-10 min)...\\n")
                    BiocManager::install(pkg, update = FALSE, ask = FALSE)
                    library(pkg, character.only = TRUE)
                    cat(pkg, ": OK\\n")
                })
            }
            cat("Step 3 complete\\n")
        ''', timeout=1200)
        self._r_print(out3[-500:] if len(out3) > 500 else out3)
        if err3:
            self._r_print(f"STDERR: {err3[-300:]}")

        # Step 4: Verify
        self._r_print("\n[Step 4] Verifying installations...")
        for pkg_name, attr in [("Seurat", "_seurat_ok"), ("monocle", "_monocle_ok")]:
            rc_v, out_v, _ = self._run_r_script(f'''
                library({pkg_name})
                cat("{pkg_name}: Ready\\n")
                cat(packageVersion("{pkg_name}"))
            ''', timeout=30)
            ok = "Ready" in out_v
            setattr(self, attr, ok)
            self._r_print(f"  {pkg_name}: {'OK' if ok else 'FAILED'}")

        # Step 5: EigenMS
        self._r_print("\n[Step 5] Installing EigenMS...")
        # Try multiple installation sources
        rc5, out5, err5 = self._run_r_script('''
            options(repos = c(CRAN = "https://cran.r-project.org"))
            library(remotes)

            tryCatch({
                library(EigenMS)
                cat("EigenMS: already installed\\n")
            }, error = function(e) {
                # Try 1: GitHub (mtmorgan)
                cat("Trying GitHub (mtmorgan/EigenMS)...\\n")
                tryCatch({
                    remotes::install_github("mtmorgan/EigenMS", upgrade = "never")
                    library(EigenMS)
                    cat("EigenMS: OK (from mtmorgan)\\n")
                }, error = function(e1) {
                    cat("mtmorgan/EigenMS failed:", conditionMessage(e1), "\\n")
                    # Try 2: BiocManager
                    cat("Trying BiocManager...\\n")
                    tryCatch({
                        library(BiocManager)
                        BiocManager::install("EigenMS", update = FALSE, ask = FALSE)
                        library(EigenMS)
                        cat("EigenMS: OK (from BiocManager)\\n")
                    }, error = function(e2) {
                        cat("BiocManager failed:", conditionMessage(e2), "\\n")
                        # Try 3: Alternative GitHub
                        cat("Trying alternative source...\\n")
                        tryCatch({
                            remotes::install_github("variani/EigenMS", upgrade = "never")
                            library(EigenMS)
                            cat("EigenMS: OK (from variani)\\n")
                        }, error = function(e3) {
                            cat("All EigenMS sources failed.\\n")
                            cat("ERROR:", conditionMessage(e3), "\\n")
                        })
                    })
                })
            })
            cat("Step 5 complete\\n")
        ''', timeout=300)
        self._r_print(out5[-500:] if len(out5) > 500 else out5)
        if err5:
            self._r_print(f"STDERR: {err5[-300:]}")
        self._eigenms_ok = "EigenMS: OK" in out5
        self._r_print(f"  EigenMS: {'OK' if self._eigenms_ok else 'FAILED — will use built-in SVD fallback'}")

        self._r_print("\n" + "=" * 50)
        self._r_print(f"Result: {self.status_message}")
        self._r_print("=" * 50)
        return self.is_ready

    # ============================================================
    # Seurat
    # ============================================================

    def run_seurat(
        self, matrix_csv: str, output_dir: str,
        n_pcs: int = 30, resolution: float = 0.8,
        min_cells: int = 3, min_features: int = 200,
    ) -> bool:
        """Run Seurat clustering + PCA + t-SNE + UMAP."""
        if not self._seurat_ok:
            self._r_print("ERROR: Seurat not installed. Run ensure_r_packages() first.")
            return False

        os.makedirs(output_dir, exist_ok=True)
        m = os.path.abspath(matrix_csv).replace("\\", "/")
        o = os.path.abspath(output_dir).replace("\\", "/")

        self._r_print(f"Running Seurat: {m} → {o}")
        rc, out, err = self._run_r_script(f'''
            library(Seurat)
            library(ggplot2)

            mat <- read.csv("{m}", row.names = 1)
            mat <- as.matrix(mat)
            cat("Loaded:", nrow(mat), "features x", ncol(mat), "cells\\n")

            sobj <- CreateSeuratObject(counts = mat, project = "SCMetabolo",
                min.cells = {min_cells}, min.features = {min_features})
            cat("Object created:", ncol(sobj), "cells x", nrow(sobj), "features\\n")

            sobj <- NormalizeData(sobj, normalization.method = "LogNormalize")
            sobj <- FindVariableFeatures(sobj, nfeatures = 2000)
            sobj <- ScaleData(sobj)
            sobj <- RunPCA(sobj, npcs = {n_pcs})
            sobj <- FindNeighbors(sobj, dims = 1:min({n_pcs}, ncol(sobj)-1))
            sobj <- FindClusters(sobj, resolution = {resolution})

            # t-SNE: use tryCatch + dynamic perplexity (Rtsne requires perplexity < n_samples/3)
            tryCatch({{
                n_cells <- ncol(sobj)
                tsne_perplexity <- min(30, max(2, floor((n_cells - 1) / 3)))
                cat("Running t-SNE with perplexity =", tsne_perplexity, "(cells =", n_cells, ")\\n")
                sobj <- RunTSNE(sobj, dims = 1:min({n_pcs}, n_cells-1),
                    perplexity = tsne_perplexity)
                cat("t-SNE done\\n")
            }}, error = function(e) {{
                cat("t-SNE skipped:", conditionMessage(e), "\\n")
            }})

            tryCatch({{
                sobj <- RunUMAP(sobj, dims = 1:min({n_pcs}, ncol(sobj)-1))
            }}, error = function(e) {{ cat("UMAP skipped\\n") }})

            setwd("{o}")
            write.csv(sobj@meta.data, "seurat_clusters.csv")
            write.csv(Embeddings(sobj, "pca"), "seurat_pca.csv")
            tryCatch({{
                write.csv(Embeddings(sobj, "tsne"), "seurat_tsne.csv")
            }}, error = function(e) {{ cat("t-SNE CSV skipped\\n") }})
            tryCatch({{
                write.csv(Embeddings(sobj, "umap"), "seurat_umap.csv")
            }}, error = function(e) {{}})

            pdf("seurat_plots.pdf", width = 10, height = 8)
            print(DimPlot(sobj, reduction = "pca", label = TRUE) + ggtitle("PCA"))
            tryCatch({{
                print(DimPlot(sobj, reduction = "tsne", label = TRUE) + ggtitle("t-SNE"))
            }}, error = function(e) {{ cat("t-SNE plot skipped\\n") }})
            tryCatch({{
                print(DimPlot(sobj, reduction = "umap", label = TRUE) + ggtitle("UMAP"))
            }}, error = function(e) {{}})
            print(ElbowPlot(sobj))
            dev.off()
            cat("Seurat done\\n")
        ''', timeout=600)

        self._r_print(out[-500:] if len(out) > 500 else out)
        if err:
            self._r_print(f"STDERR: {err[-300:]}")
        return "Seurat done" in out

    # ============================================================
    # Monocle
    # ============================================================

    def run_monocle(
        self, matrix_csv: str, output_dir: str, max_components: int = 2
    ) -> bool:
        """Run Monocle pseudotime analysis."""
        if not self._monocle_ok:
            self._r_print("ERROR: Monocle not installed. Run ensure_r_packages() first.")
            return False

        os.makedirs(output_dir, exist_ok=True)
        m = os.path.abspath(matrix_csv).replace("\\", "/")
        o = os.path.abspath(output_dir).replace("\\", "/")

        self._r_print(f"Running Monocle: {m} → {o}")
        rc, out, err = self._run_r_script(f'''
            library(monocle)

            mat <- read.csv("{m}", row.names = 1)
            mat <- as.matrix(mat)
            cat("Loaded:", nrow(mat), "x", ncol(mat), "\\n")

            pd <- new("AnnotatedDataFrame", data = data.frame(
                row.names = colnames(mat), sample = colnames(mat)))
            fd <- new("AnnotatedDataFrame", data = data.frame(
                row.names = rownames(mat), gene_short_name = rownames(mat)))

            # Use gaussianff() for continuous metabolomics data
            cds <- newCellDataSet(as(mat, "sparseMatrix"),
                phenoData = pd, featureData = fd,
                expressionFamily = gaussianff())
            cds <- estimateSizeFactors(cds)

            # Select ordering genes by variance (not dispersion)
            feature_var <- apply(exprs(cds), 1, var, na.rm = TRUE)
            n_ordering <- min(500, length(feature_var))
            og <- names(sort(feature_var, decreasing = TRUE))[1:n_ordering]
            cat("Ordering genes selected:", length(og), "/", length(feature_var), "\\n")
            cds <- setOrderingFilter(cds, og)

            cds <- reduceDimension(cds, max_components = {max_components},
                reduction_method = "DDRTree", norm_method = "none")

            cat("MANUAL-PSEUDOTIME-V3\\n")
            cat("Building MST from reducedDim coordinates...\\n")

            # DDRTree output structure changed in newer DDRTree versions:
            # pr_graph_cell_proj_tree / Z no longer exist in auxOrderingData.
            # We build our own MST from reducedDimS(cds) coordinates.
            coords <- t(reducedDimS(cds))  # n_cells x n_dim
            cell_names <- colnames(cds)

            # Compute distance matrix and build MST using new igraph API
            dmat <- as.matrix(dist(coords))
            g <- igraph::graph_from_adjacency_matrix(dmat, mode = "undirected",
                weighted = TRUE, diag = FALSE)
            mst_tree <- igraph::mst(g)  # new igraph API (not minimum.spanning.tree)

            # Root cell: closest to centroid
            centroid <- colMeans(coords)
            dist_to_centroid <- rowSums((coords - matrix(centroid,
                nrow = nrow(coords), ncol = ncol(coords), byrow = TRUE))^2)
            root_idx <- which.min(dist_to_centroid)
            cat("Root cell:", cell_names[root_idx], "\\n")

            # Pseudotime = shortest-path distance from root in MST
            all_dist <- as.numeric(igraph::distances(mst_tree, v = root_idx))
            pseudotime <- all_dist
            cat("Pseudotime range:", min(pseudotime), "-", max(pseudotime), "\\n")

            pData(cds)$Pseudotime <- pseudotime
            n_states <- min(3, max(1, length(unique(round(pseudotime, 4)))))
            pData(cds)$State <- as.numeric(cut(pseudotime,
                breaks = n_states, labels = FALSE, include.lowest = TRUE))
            pData(cds)$State[is.na(pData(cds)$State)] <- 1
            cat("States:", paste(unique(pData(cds)$State), collapse=","), "\\n")

            setwd("{o}")
            write.csv(pData(cds), "monocle_pseudotime.csv")
            write.csv(coords, "monocle_coords.csv")
            write.csv(og, "monocle_ordering_genes.csv")

            # Custom trajectory plot — plot_cell_trajectory() is broken
            # (uses dplyr::select_() which is defunct in dplyr >= 0.7.0)
            library(ggplot2)
            df <- data.frame(
                Dim1 = coords[, 1], Dim2 = coords[, 2],
                Pseudotime = pData(cds)$Pseudotime,
                State = factor(pData(cds)$State))
            edges <- igraph::as_edgelist(mst_tree)
            edge_df <- data.frame(
                x = coords[edges[, 1], 1], y = coords[edges[, 1], 2],
                xend = coords[edges[, 2], 1], yend = coords[edges[, 2], 2])

            p1 <- ggplot(df, aes(x = Dim1, y = Dim2)) +
                geom_segment(data = edge_df, aes(x = x, y = y,
                    xend = xend, yend = yend),
                    color = "grey70", linewidth = 0.3) +
                geom_point(aes(color = Pseudotime), size = 1.5) +
                scale_color_gradient(low = "lightblue", high = "darkred") +
                labs(title = "Monocle Trajectory — Pseudotime",
                    x = "Component 1", y = "Component 2") +
                theme_minimal()

            p2 <- ggplot(df, aes(x = Dim1, y = Dim2)) +
                geom_segment(data = edge_df, aes(x = x, y = y,
                    xend = xend, yend = yend),
                    color = "grey70", linewidth = 0.3) +
                geom_point(aes(color = State), size = 1.5) +
                labs(title = "Monocle Trajectory — State",
                    x = "Component 1", y = "Component 2") +
                theme_minimal()

            pdf("monocle_plots.pdf", width = 10, height = 8)
            print(p1)
            print(p2)
            dev.off()
            cat("Monocle done\\n")
        ''', timeout=600)

        self._r_print(out[-500:] if len(out) > 500 else out)
        if err:
            self._r_print(f"STDERR: {err[-300:]}")
        return "Monocle done" in out

    # ============================================================
    # EigenMS
    # ============================================================

    def run_eigenms(
        self, matrix_csv: str, output_dir: str, n_eig: int = 10
    ) -> bool:
        """Run EigenMS SVD normalization (with fallback to base SVD)."""
        os.makedirs(output_dir, exist_ok=True)
        m = os.path.abspath(matrix_csv).replace("\\", "/")
        o = os.path.abspath(output_dir).replace("\\", "/")

        if self._eigenms_ok:
            self._r_print(f"Running EigenMS: {m} → {o}")
            rc, out, err = self._run_r_script(f'''
                library(EigenMS)
                mat <- read.csv("{m}", row.names = 1)
                mat <- as.matrix(mat)
                mat[mat == 0] <- NA
                mat <- t(mat)
                mat <- log2(mat)

                result <- eig_norm1(m = mat, treatment = NULL,
                    prot.info = NULL, write_to_file = "")

                norm <- 2^result$norm_matrix
                setwd("{o}")
                write.csv(t(norm), "eigenms_normalized.csv")
                write.csv(result$eigentrends, "eigenms_eigentrends.csv")

                pdf("eigenms_diagnostics.pdf", width = 10, height = 6)
                plot(result$eigentrends, type = "b",
                    main = "Eigentrends", xlab = "Eigenvalue", ylab = "Trend")
                dev.off()
                cat("EigenMS done\\n")
            ''', timeout=300)
        else:
            # Fallback: base SVD
            self._r_print(f"Running SVD normalization (EigenMS not installed): {m} → {o}")
            rc, out, err = self._run_r_script(f'''
                mat <- read.csv("{m}", row.names = 1)
                mat <- as.matrix(mat)
                mat[mat == 0] <- NA

                log_mat <- log2(mat)
                for (i in 1:ncol(log_mat)) {{
                    log_mat[is.na(log_mat[, i]), i] <- mean(log_mat[, i], na.rm = TRUE)
                }}

                svd_res <- svd(scale(log_mat, scale = FALSE))
                k <- min(3, ncol(svd_res$v))
                corrected <- log_mat - svd_res$u[, 1:k] %*%
                    diag(svd_res$d[1:k]) %*% t(svd_res$v[, 1:k])
                normalized <- 2^corrected

                setwd("{o}")
                write.csv(normalized, "eigenms_normalized.csv")

                pdf("eigenms_diagnostics.pdf", width = 10, height = 6)
                plot(svd_res$d, type = "b", main = "Singular Values",
                    xlab = "Component", ylab = "Value")
                dev.off()
                cat("SVD done\\n")
            ''', timeout=300)

        self._r_print(out[-500:] if len(out) > 500 else out)
        if err:
            self._r_print(f"STDERR: {err[-300:]}")
        return "EigenMS done" in out or "SVD done" in out


# ============================================================
# Singleton
# ============================================================

_bridge_instance: Optional[RBridge] = None


def get_r_bridge() -> RBridge:
    """Get or create the singleton RBridge instance."""
    global _bridge_instance
    if _bridge_instance is None:
        _bridge_instance = RBridge()
    return _bridge_instance
