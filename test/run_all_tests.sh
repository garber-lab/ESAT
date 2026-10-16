#!/bin/bash
# Run all ESAT tests against a built jar (default: newest dist/esat.*.jar).
# Requirements: java (>=8), samtools, python3. Each test creates and cleans up its own temp directory.
set -euo pipefail
cd "$(dirname "$0")/.."
JAR=${1:-$(ls -t dist/esat.*.jar | head -1)}
echo "Testing $JAR"
echo "== 1. Golden regression vs production umihack jar (default and -umiMode gene)"
python3 test/regression/test_regression.py "$JAR"
echo "== 2. Sensitivity: -umiMode position must change the outputs"
python3 test/regression/test_regression.py "$JAR" --expect-fail -umiMode position
echo "== 3. Bulk UMI feature tests (-umiMode / -umiDropN)"
python3 test/bulkumi/test_bulkumi.py "$JAR"
echo "== 4. UMI saturation correction (-umiSaturation)"
python3 test/bulkumi/test_saturation.py "$JAR"
echo "ALL ESAT TESTS PASSED"
