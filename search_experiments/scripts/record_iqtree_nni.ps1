<#
    Records IQ-TREE's own built-in stochastic NNI tree search (the real
    pipeline -- main/phyloanalysis.cpp + tree/iqtree.cpp -- run via the full
    iqtree3.exe, NOT spr_topology_test's standalone SPR API) into a CSV with
    the same columns spr_topology_test's --hillclimb "record" flag writes:
    run_id,candidates,time_elapsed,logL,true_minus_current
    (see recordSpreadsheetPath/appendRecordRow in tree/spr_topology_test.cpp).

    Two things distinguish this from just eyeballing an `iqtree3 -s ... -m
    ...` log by hand:

    1) "Own start": the search NEVER receives the AliSim ground-truth tree
       as a starting tree or constraint (no -t/-g) -- it builds its own
       BIONJ/parsimony start, exactly like a real analysis would. The true
       tree is only ever used in a SEPARATE `-te` (fixed-tree) evaluation
       run, purely to compute the true_minus_current reference column, the
       same "known answer, only used to score the attempt" role
       sim.treefile plays in spr_topology_test's own --hillclimb.

    2) Full wall-clock timing. IQ-TREE prints "Iteration N / LogL: ... /
       Time: HhMmSs" per iteration and "Total wall-clock time used: ..." at
       the end, but BOTH are measured from params.start_real_time, which
       phyloanalysis.cpp only sets AFTER alignment reading, PLL setup, and
       the initial BIONJ/parsimony distance computation have already run
       (see the runPhyloAnalysis call order around
       "params.start_real_time = getRealTime();") -- i.e. IQ-TREE's own
       clock starts only once a start tree already exists. This script
       instead times the whole external process (spawn to exit) and adds
       the (externalTotal - internalTotal) difference as a constant offset
       to every per-iteration timestamp, so time_elapsed in the CSV really
       does span from process launch to completion.

    Usage (mirrors spr_topology_test --hillclimb's calling convention):
        record_iqtree_nni.ps1 <alisim-tree.treefile> <nstop> [gtr] [quiet]

    <alisim-tree.treefile>  path to the AliSim ground-truth tree; the
                             alignment is derived by replacing ".treefile"
                             with ".fa", same convention as --hillclimb.
    <nstop>                 passed straight through as iqtree3's own -nstop
                             (number of unsuccessful iterations before the
                             stochastic NNI search gives up) -- the closest
                             built-in-NNI analogue of --hillclimb's
                             <max-steps>.
    gtr                     search under GTR+FO (ML-estimated rates/freqs)
                             instead of plain JC.
    quiet                   only print the run summary, not one line per
                             recorded iteration (does not affect what's
                             written to the CSV).

    Example:
        test_scripts/record_iqtree_nni.ps1 sim.treefile 100
        test_scripts/record_iqtree_nni.ps1 sim.treefile 100 gtr quiet
#>

param(
    [Parameter(Mandatory=$true, Position=0)][string]$TreeFile,
    [Parameter(Mandatory=$true, Position=1)][int]$MaxSteps,
    [Parameter(Position=2, ValueFromRemainingArguments=$true)][string[]]$Flags,
    [string]$IQTreeBin = "build/iqtree3.exe",
    [string]$WorkDir = (Join-Path $env:TEMP "iqtree_nni_runs"),
    [string]$CsvPathOverride = ""
)

$ErrorActionPreference = "Stop"

$UseGtr = $Flags -contains "gtr"
$Quiet  = $Flags -contains "quiet"

if (-not (Test-Path $TreeFile)) {
    throw "Tree file not found: $TreeFile"
}
if (-not (Test-Path $IQTreeBin)) {
    throw "iqtree3 binary not found at $IQTreeBin -- build it first: cmake --build build --target iqtree3"
}

$Alignment = $TreeFile -replace '\.treefile$', '.fa'
if ($Alignment -eq $TreeFile) {
    throw "Expected a .treefile path (alignment is derived by replacing .treefile with .fa): $TreeFile"
}
if (-not (Test-Path $Alignment)) {
    throw "Alignment not found next to tree file (expected $Alignment)"
}

$Model = if ($UseGtr) { "GTR+FO" } else { "JC" }

New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null

$Timestamp = Get-Date -Format "yyyyMMdd-HHmmss"
$RunId     = "${Timestamp}_iqtree_nni_ownstart_nstop${MaxSteps}"
$Prefix    = Join-Path $WorkDir $RunId

function Add-CsvLine {
    param([string]$Path, [string]$Line)
    [System.IO.File]::AppendAllText($Path, "$Line`r`n")
}

# --- Step 1: evaluate the true AliSim tree's own logL under the same model.
# This run is NEVER fed into the search below; it exists only to compute
# the true_minus_current reference column, same role sim.treefile plays in
# spr_topology_test's --hillclimb.
# iqtree3 tees everything it prints into its own <prefix>.log file
# regardless of what its stdout is connected to, so unless quiet was
# requested we don't redirect stdout at all -- with -NoNewWindow that lets
# it inherit this console directly and stream live, instead of the whole
# run going silent until it exits. Post-run parsing below always reads the
# authoritative <prefix>.log file, never this stdout stream.
function Invoke-IQTree {
    param([string[]]$IqArgs)
    if ($Quiet) {
        return Start-Process -FilePath $IQTreeBin -ArgumentList $IqArgs -NoNewWindow -PassThru -Wait -RedirectStandardOutput "NUL"
    } else {
        return Start-Process -FilePath $IQTreeBin -ArgumentList $IqArgs -NoNewWindow -PassThru -Wait
    }
}

if (-not $Quiet) { Write-Host "Evaluating true tree's likelihood (reference only -- not used by the search)..." }

$trueEvalPrefix = "${Prefix}_trueeval"
$trueEvalArgs = @("-s", $Alignment, "-te", $TreeFile, "-m", $Model, "--show-lh", "-T", "1", "--prefix", $trueEvalPrefix, "--redo")
$p = Invoke-IQTree -IqArgs $trueEvalArgs
if ($p.ExitCode -ne 0) {
    throw "trueeval run failed (exit $($p.ExitCode)); see $trueEvalPrefix.log"
}

$trueIqtreeFile = "$trueEvalPrefix.iqtree"
$trueMatch = Select-String -Path $trueIqtreeFile -Pattern "Log-likelihood of the tree:\s*(-?[0-9.]+)" | Select-Object -First 1
if (-not $trueMatch) {
    throw "Could not find true tree log-likelihood in $trueIqtreeFile"
}
$TrueLogL = [double]$trueMatch.Matches[0].Groups[1].Value
if (-not $Quiet) { Write-Host "True tree logL: $TrueLogL" }

# --- Step 2: the actual search, own BIONJ/parsimony start (no -t/-g), timed
# externally from just before the process is spawned to just after it exits.
if (-not $Quiet) { Write-Host "Running iqtree3's own NNI search (own start, model=$Model, nstop=$MaxSteps)..." }

$mainPrefix = "${Prefix}_run"
$mainArgs = @("-s", $Alignment, "-m", $Model, "-v", "-T", "1", "--prefix", $mainPrefix, "--redo", "-nstop", "$MaxSteps")

$externalStart = Get-Date
$p = Invoke-IQTree -IqArgs $mainArgs
$externalEnd = Get-Date
if ($p.ExitCode -ne 0) {
    throw "NNI search run failed (exit $($p.ExitCode)); see $mainPrefix.log"
}
$externalTotalSec = ($externalEnd - $externalStart).TotalSeconds

# --- Step 3: parse the log for per-iteration progress and IQ-TREE's own
# (partial) total, then correct every timestamp by the difference.
$logFile = "$mainPrefix.log"
$logLines = Get-Content $logFile

$totalMatch = $logLines | Select-String "Total wall-clock time used:\s*([0-9.]+)\s*sec" | Select-Object -Last 1
if (-not $totalMatch) {
    throw "Could not find 'Total wall-clock time used' in $logFile"
}
$InternalTotalSec = [double]$totalMatch.Matches[0].Groups[1].Value
$Offset = $externalTotalSec - $InternalTotalSec

if (-not $Quiet) {
    Write-Host ("External process time: {0:N2}s -- IQ-TREE-reported time: {1:N2}s -- pre-search-timer overhead added back: {2:N2}s" `
        -f $externalTotalSec, $InternalTotalSec, $Offset)
}

$iterPattern = '^(?:Iteration|Bootstrap) (\d+) / LogL: (-?[0-9.]+) / Time: (\d+)h:(\d+)m:(\d+)s'
$rows = @()
foreach ($line in $logLines) {
    $m = [regex]::Match($line, $iterPattern)
    if ($m.Success) {
        $candidates  = [long]$m.Groups[1].Value
        $logL        = [double]$m.Groups[2].Value
        $internalSec = [double]$m.Groups[3].Value * 3600 + [double]$m.Groups[4].Value * 60 + [double]$m.Groups[5].Value
        $rows += [PSCustomObject]@{
            Candidates  = $candidates
            TimeElapsed = [math]::Round($internalSec + $Offset, 3)
            LogL        = $logL
        }
    }
}
if ($rows.Count -eq 0) {
    throw "No 'Iteration N / LogL: ... / Time: ...' lines found in $logFile -- was -v honored?"
}

# --- Step 4: append to the record spreadsheet, same format as
# recordSpreadsheetPath/appendRecordRow in tree/spr_topology_test.cpp.
$sanitizedModel = ($Model -replace '[^A-Za-z0-9]', '_')
$csvPath = if ($CsvPathOverride) { $CsvPathOverride } else { "record_${sanitizedModel}_iqtree_nni_ownstart.csv" }
$needsHeader = (-not (Test-Path $csvPath)) -or ((Get-Item $csvPath).Length -eq 0)
if ($needsHeader) {
    Add-CsvLine $csvPath "run_id,candidates,time_elapsed,logL,true_minus_current"
}

foreach ($row in $rows) {
    $trueMinusCurrent = $TrueLogL - $row.LogL
    Add-CsvLine $csvPath "$RunId,$($row.Candidates),$($row.TimeElapsed),$($row.LogL),$trueMinusCurrent"
    if (-not $Quiet) {
        Write-Host ("  candidates={0} time_elapsed={1}s logL={2} true_minus_current={3}" `
            -f $row.Candidates, $row.TimeElapsed, $row.LogL, $trueMinusCurrent)
    }
}

Write-Host "Appended $($rows.Count) rows to $csvPath (run_id=$RunId)"
Write-Host "Raw iqtree3 outputs kept under $WorkDir (prefix $RunId*) for inspection."
