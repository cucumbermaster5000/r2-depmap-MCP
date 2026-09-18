$ErrorActionPreference = 'Stop'
$rCommand = Get-Command Rscript.exe -ErrorAction SilentlyContinue
if ($rCommand) {
    $rExecutable = $rCommand.Source
} else {
    $rExecutable = Get-ChildItem -Path "$env:ProgramFiles\R\R-*\bin\Rscript.exe" -ErrorAction SilentlyContinue |
        Sort-Object { [version]($_.Directory.Parent.Name -replace '^R-', '') } -Descending |
        Select-Object -First 1 -ExpandProperty FullName
}
if (-not $rExecutable) { throw 'Install R or add Rscript.exe to PATH.' }
& $rExecutable --vanilla (Join-Path $PSScriptRoot '..\R\server.R')
exit $LASTEXITCODE
