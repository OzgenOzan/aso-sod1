"""
Master Pipeline Runner
========================
Executes all phases sequentially.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def run_pipeline():
    start_time = time.time()
    print("=" * 70)
    print("  SOD1 ASO INHIBITION PREDICTION PIPELINE")
    print("  Master Pipeline Runner")
    print("=" * 70)
    print()

    # Phase 1
    from phase1_audit import main as phase1
    phase1()

    # Phase 2
    from phase2_cleaning import main as phase2
    phase2()

    # Phase 3
    from phase3_clustering import main as phase3
    phase3()

    # Phase 4
    from phase4_features import main as phase4
    phase4()

    # Phase 5
    from phase5_split import main as phase5
    phase5()

    # Phase 6
    from phase6_external import main as phase6
    phase6()

    # Phase 7
    from phase7_modeling import main as phase7
    phase7()

    # Phase 8
    from phase8_validation import main as phase8
    phase8()

    # Phase 9
    from phase9_tofersen import main as phase9
    phase9()

    elapsed = time.time() - start_time
    print("=" * 70)
    print(f"  ALL PHASES COMPLETE -- Total time: {elapsed:.1f}s")
    print("=" * 70)


if __name__ == "__main__":
    run_pipeline()
