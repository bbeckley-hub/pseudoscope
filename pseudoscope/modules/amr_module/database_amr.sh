#!/bin/bash
# database_amr.sh - Fetch AMRFinder binaries and download the latest database
set -euo pipefail

REPO_URL="https://github.com/bbeckley-hub/amr.git"
BRANCH="master"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="$SCRIPT_DIR/data/amrfinder_db"

NCBI_MIRROR="${AMR_NCBI_MIRROR:-https://ftp.ncbi.nlm.nih.gov}"
OFFLINE="${AMR_OFFLINE:-0}"
MAX_ATTEMPTS=3

log()  { echo "[$(date +%H:%M:%S)] $*"; }
warn() { echo "[$(date +%H:%M:%S)] WARN: $*" >&2; }
die()  { echo "[$(date +%H:%M:%S)] ERROR: $*" >&2; exit 1; }

log "Fetching AMR binaries and helper tools from $REPO_URL"
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

git clone --depth 1 --branch "$BRANCH" "$REPO_URL" "$TEMP_DIR"

for dir in bin stx; do
    if [ -d "$TEMP_DIR/$dir" ]; then
        log "Updating $dir/"
        rm -rf "$SCRIPT_DIR/$dir"
        cp -r "$TEMP_DIR/$dir" "$SCRIPT_DIR/"
    fi
done

chmod -R a+x "$SCRIPT_DIR/bin/" 2>/dev/null || true
mkdir -p "$DB_DIR"

[ -x "$SCRIPT_DIR/bin/amrfinder_update" ] \
    || die "amrfinder_update binary not found at $SCRIPT_DIR/bin/"

if [ "$OFFLINE" = "1" ]; then
    ls -d "$DB_DIR"/20*/ >/dev/null 2>&1 \
        || die "Offline mode requested but no existing database found in $DB_DIR"
    LATEST=$(ls -d "$DB_DIR"/20*/ | sort | tail -1)
    log "Offline mode: using existing database $LATEST"
    exit 0
fi

PROBE_URL="$NCBI_MIRROR/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/"
log "Probing NCBI endpoint: $PROBE_URL"
if ! curl -sSfL --connect-timeout 10 --max-time 20 -o /dev/null "$PROBE_URL"; then
    warn "Cannot reach NCBI at $NCBI_MIRROR"
    warn ""
    warn "Diagnostic steps (try each in order):"
    warn "  1. Direct probe:        curl -I $NCBI_MIRROR/"
    warn "  2. Force IPv4:          curl -4 -I $NCBI_MIRROR/"
    warn "  3. Via proxy:           export HTTPS_PROXY=http://your-proxy:port"
    warn "  4. Check NCBI status:   https://www.ncbi.nlm.nih.gov/"
    warn "  5. Check DNS:           nslookup ftp.ncbi.nlm.nih.gov"
    warn ""
    warn "Workarounds:"
    warn "  a. Install the bioconda database snapshot (no NCBI access needed):"
    warn "       micromamba install -n pseudoscope -c bioconda amrfinderplus-db"
    warn "     Then copy or symlink its data into:"
    warn "       $DB_DIR"
    warn ""
    warn "  b. Use an existing local database (already downloaded):"
    warn "       AMR_OFFLINE=1 bash database_amr.sh"
    warn ""
    warn "  c. Point at a mirror:"
    warn "       AMR_NCBI_MIRROR=https://your-mirror.example.org bash database_amr.sh"
    die "NCBI unreachable"
fi

log "Downloading the latest AMR database (this may take a few minutes)..."
for attempt in $(seq 1 "$MAX_ATTEMPTS"); do
    log "Attempt $attempt/$MAX_ATTEMPTS"
    if "$SCRIPT_DIR/bin/amrfinder_update" --database "$DB_DIR"; then
        break
    fi
    if [ "$attempt" -lt "$MAX_ATTEMPTS" ]; then
        backoff=$((attempt * 15))
        warn "Attempt $attempt failed, retrying in ${backoff}s..."
        sleep "$backoff"
    else
        die "All $MAX_ATTEMPTS attempts to download AMR database failed"
    fi
done

ls -d "$DB_DIR"/20*/ >/dev/null 2>&1 \
    || die "amrfinder_update returned success but no database folder (20*) exists in $DB_DIR"

LATEST=$(ls -d "$DB_DIR"/20*/ | sort | tail -1)
SIZE=$(du -sh "$LATEST" | cut -f1)
VERSION=$(cat "$LATEST/version.txt" 2>/dev/null || basename "$LATEST")
log "Database installed: $LATEST"
log "  Version: $VERSION"
log "  Size:    $SIZE"
log "AMR setup complete."
