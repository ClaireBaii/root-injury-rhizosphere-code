# Check if script is running in PowerShell
if ($PSVersionTable.PSVersion.Major -lt 5) {
  Write-Error "This script requires PowerShell 5.0 or later."
  exit 1
}

$ErrorActionPreference = "Stop"

# Define environment name
$ENV_NAME = "paper_fig_env"

# Use absolute path to avoid PATH issues (as seen in legacy code)
$CONDA_EXE = "D:\conda\condabin\conda.bat"
if (-not (Test-Path $CONDA_EXE)) {
  # Fallback to 'conda' if specific path doesn't exist
  $CONDA_EXE = "conda" 
}

# 1. Create environment if it doesn't exist (update if it does)
Write-Host "Ensuring conda environment '$ENV_NAME' is installed..." -ForegroundColor Cyan

& $CONDA_EXE env update -f environment.yaml --prune
if ($LASTEXITCODE -ne 0) {
  Write-Error "Conda environment update failed."
  exit $LASTEXITCODE
}

# 2. Verify installation using 'conda run'
Write-Host "----------------------------------------------------------------" -ForegroundColor Green
Write-Host "Verifying R packages in '$ENV_NAME'..." -ForegroundColor Green
Write-Host "----------------------------------------------------------------" -ForegroundColor Green

$r_script = @"
pkgs <- c(
  'ggplot2', 'igraph', 'Hmisc', 'ggraph', 'RColorBrewer', 'scales',
  'circlize', 'vegan', 'data.table', 'mixOmics', 'ComplexHeatmap',
  'pheatmap', 'ggrepel'
)

message('Checking package status...')
installed <- sapply(pkgs, requireNamespace, quietly = TRUE)
if (all(installed)) {
  message('[OK] All required packages are found!')
} else {
  missing <- pkgs[!installed]
  stop('Missing packages: ', paste(missing, collapse = ', '))
}
"@

# Verify using user-requested temporary file logic
$tmp_r = Join-Path $env:TEMP "verify_packages_${ENV_NAME}.R"
Set-Content -Path $tmp_r -Value $r_script -Encoding UTF8

try {
  Write-Host "Running verification script..." -ForegroundColor Gray
    
  # Run using the configured conda executable
  & $CONDA_EXE run -n $ENV_NAME Rscript $tmp_r
    
  # Capture exit code
  $VerifyExitCode = $LASTEXITCODE
}
finally {
  # Clean up
  Remove-Item $tmp_r -ErrorAction SilentlyContinue
}

# 3. Final instructions
if ($VerifyExitCode -eq 0) {
  Write-Host "----------------------------------------------------------------" -ForegroundColor Cyan
  Write-Host "Setup Successful!" -ForegroundColor Cyan
  Write-Host "To start using the environment, run:" -ForegroundColor White
  Write-Host "    conda activate $ENV_NAME" -ForegroundColor Yellow
  Write-Host "----------------------------------------------------------------" -ForegroundColor Cyan
}
else {
  Write-Host "----------------------------------------------------------------" -ForegroundColor Red
  Write-Warning "Verification Failed. Checks above for errors."
  Write-Host "----------------------------------------------------------------" -ForegroundColor Red
  exit $VerifyExitCode
}
$ErrorActionPreference = "Stop"

# Define environment name
$ENV_NAME = "paper_fig_env"

# 1. Create environment if it doesn't exist (update if it does)
# 使用绝对路径以避免 PATH 变量未生效的问题
$CONDA_EXE = "D:\conda\condabin\conda.bat"

Write-Host "Ensuring conda environment '$ENV_NAME' is installed..."
# 使用 cmd /c 运行 conda 以确保在 PowerShell 中也能正确捕获批处理文件的执行
& $CONDA_EXE env update -f environment.yaml --prune
if ($LASTEXITCODE -ne 0) {
  Write-Error "Conda environment update failed."
}

# 2. Verify installation using 'conda run'
Write-Host "----------------------------------------------------------------"
Write-Host "Verifying R packages in '$ENV_NAME'..."
Write-Host "----------------------------------------------------------------"

$r_script = @"
pkgs <- c(
  'ggplot2', 'igraph', 'Hmisc', 'ggraph', 'RColorBrewer', 'scales',
  'circlize', 'vegan', 'data.table', 'mixOmics', 'ComplexHeatmap',
  'pheatmap', 'ggrepel'
)

message('Checking package status...')
installed <- sapply(pkgs, requireNamespace, quietly = TRUE)
if (all(installed)) {
  message('[OK] All required packages are found!')
} else {
  missing <- pkgs[!installed]
  stop('Missing packages: ', paste(missing, collapse = ', '))
}
"@

# 将 R 代码保存为临时文件以避免引号转义问题，或者直接传递
# 在 PowerShell 中传递多行字符串给 cmd/conda 可能很棘手，这里我们尝试尽量简洁
# 注意：conda run 在 Windows 上有时行为怪异，我们构建单行命令
$r_one_line = $r_script -replace "`r`n", ";" -replace "`n", ";" -replace '"', '\"'

& $CONDA_EXE run -n $ENV_NAME Rscript -e "$r_one_line"

# 3. Final instructions
if ($LASTEXITCODE -eq 0) {
  Write-Host "----------------------------------------------------------------"
  Write-Host "Setup Successful!"
  Write-Host "To start using the environment, run:"
  Write-Host "    conda activate $ENV_NAME"
  Write-Host "----------------------------------------------------------------"
}
else {
  Write-Host "----------------------------------------------------------------"
  Write-Warning "Verification Failed. Checks above for errors."
  Write-Host "----------------------------------------------------------------"
}
