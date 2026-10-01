"""Flight data CSV reader for fin_flutter analysis pipeline."""

import csv
from pathlib import Path


def read_flight_data(csv_path: str | Path) -> tuple[float, float]:
    """Extract maximum velocity and the altitude at which it occurs from a flight data CSV.

    The CSV format expected is:
        # comment lines (skipped)
        Time (s), Altitude (ft), Vertical velocity (m/s)
        <float>, <float>, <float>
        ...

    Args:
        csv_path: Path to the flight data CSV file.

    Returns:
        Tuple of (max_velocity_m_s, altitude_at_max_velocity_ft) where:
            max_velocity_m_s: Maximum vertical velocity in m/s.
            altitude_at_max_velocity_ft: Altitude in ft at the time of maximum velocity.

    Raises:
        FileNotFoundError: If the CSV file does not exist.
        ValueError: If the file contains no valid data rows.
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Flight data file not found: {csv_path}")

    max_velocity: float = float("-inf")
    altitude_at_max: float = 0.0

    with csv_path.open(newline="") as f:
        reader = csv.reader(f)
        rows_read = 0
        for row in reader:
            # Skip comment lines and blank lines
            if not row or row[0].strip().startswith("#"):
                continue
            try:
                altitude = float(row[1])
                velocity = float(row[2])
            except (ValueError, IndexError):
                continue
            rows_read += 1
            if velocity > max_velocity:
                max_velocity = velocity
                altitude_at_max = altitude

    if rows_read == 0:
        raise ValueError(f"No valid data rows found in {csv_path}")

    return max_velocity, altitude_at_max
