# --- portable roots -----------------------------------------------------------------------------
# Injected by scripts/build_repo.py; see the note in the Python header. Set CKM_P4_BASE to override.
P4_BASE <- Sys.getenv("CKM_P4_BASE")
if (!nzchar(P4_BASE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .self <- sub("^--file=", "", .args[grep("^--file=", .args)])
  .here <- if (length(.self)) dirname(.self[1]) else "."
  P4_BASE <- normalizePath(file.path(.here, ".."), winslash = "/", mustWork = FALSE)
}
P4_ROOT <- Sys.getenv("CKM_P4_ROOT"); if (!nzchar(P4_ROOT)) P4_ROOT <- dirname(P4_BASE)
CKM_ROOT <- Sys.getenv("CKM_ROOT");   if (!nzchar(CKM_ROOT)) CKM_ROOT <- dirname(P4_ROOT)
SHARED_LIB <- Sys.getenv("CKM_SHARED_LIB")
if (!nzchar(SHARED_LIB)) SHARED_LIB <- file.path(CKM_ROOT, "_shared")
# ------------------------------------------------------------------------------------------------
# =============================================================================
# build_renv_lock.R — record the R environment the analysis actually ran in (round-2 item B-2 #22).
#
# Deliberately NOT `renv::init()`. Initialising renv in this project would rewrite .Rprofile and
# start relocating the user's library, which is a destructive side effect for a reproducibility
# record. Instead the lockfile is built with renv's own dependency scanner and snapshot machinery
# pointed at a throwaway project directory, so the user's library and profile are untouched.
#
# The recorded set is DISCOVERED from the scripts (renv::dependencies), not typed, so a package the
# analysis uses cannot be missing from the record because somebody forgot to add it.
# =============================================================================
BASE <- P4_BASE
OUT  <- file.path(BASE, "manifest", "renv.lock")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

stopifnot(requireNamespace("renv", quietly = TRUE))

deps <- renv::dependencies(path = c(file.path(BASE, "scripts"), file.path(BASE, "figures")),
                           errors = "ignored", progress = FALSE)
pkgs <- sort(unique(deps$Package))
pkgs <- pkgs[vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
cat(sprintf("packages discovered in the analysis scripts and installed here: %d\n", length(pkgs)))

tmp <- file.path(tempdir(), "ckm_p4_renv")
unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
renv::snapshot(project = tmp, packages = pkgs, lockfile = OUT,
               type = "packages", prompt = FALSE, force = TRUE)

lock <- jsonlite::fromJSON(OUT)
cat(sprintf("wrote %s\n", OUT))
cat(sprintf("R version recorded : %s\n", lock$R$Version))
cat(sprintf("packages recorded  : %d\n", length(lock$Packages)))

# The record is only useful if it covers what the figures and the MR battery actually load.
must <- c("ggplot2", "patchwork", "scales", "dplyr", "TwoSampleMR", "MRPRESSO")
missing <- setdiff(intersect(must, pkgs), names(lock$Packages))
if (length(missing)) {
  stop("renv.lock is missing packages the analysis uses: ", paste(missing, collapse = ", "))
}
cat("core analysis packages present in the lockfile:",
    paste(intersect(must, names(lock$Packages)), collapse = ", "), "\n")
