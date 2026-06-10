from __future__ import annotations

import argparse
from pathlib import Path

from GTsRunner import GTsRunner

def main() -> None:
    parser = argparse.ArgumentParser(description="Predict and remove solvent by GTsR.")
    parser.add_argument("--cif", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.5)
    args = parser.parse_args()

    runner = GTsRunner(checkpoint=args.checkpoint)
    result = runner.predict(
        cif=args.cif,
        output=args.output,
        threshold=args.threshold,
    )
    print(
        f"Predicted {result['num_solvent_atoms']} solvent atoms and "
        f"{result['num_framework_atoms']} framework atoms."
    )
    print(f"Framework CIF: {result['framework']}")
    if result["solvent"]:
        print(f"Solvent CIF: {result['solvent']}")

if __name__ == "__main__":
    main()
