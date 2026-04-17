from __future__ import annotations

import pandas as pd

from gatk_sv_ploidy.eval import calculate_metrics


def test_calculate_metrics_writes_report(tmp_path) -> None:
    pred_df = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3"],
            "true_aneuploidy_type": ["NORMAL", "TURNER", "NORMAL"],
            "predicted_aneuploidy_type": ["NORMAL", "NORMAL", "NORMAL"],
        }
    )

    calculate_metrics(pred_df, str(tmp_path))

    report = (tmp_path / "metrics_report.txt").read_text()
    assert "Overall Accuracy: 2/3" in report
    assert "Confusion Matrix" in report
