"""
Configuration file containing run counts per subject and session
"""

# Database of runs per session for each subject
RUNS_DATABASE = {
    'sub-001': {
        'ses-01': 16,
        'ses-02': 16,
        'ses-03': 16,
        'ses-04': 16,
        'ses-05': 5,
        'ses-06': 5,
        'ses-07': 5
    },
    'sub-002': {
        'ses-01': 9,
        'ses-02': 5
    },
    'sub-003': {
        'ses-01': 16
    },
    'sub-004': {
        'ses-01': 5,
        'ses-02': 12
    },
    'sub-005': {
        'ses-01': 11,
        'ses-02': 11,
        'ses-03': 11
    },
    'sub-006': {
        'ses-01': 5,
        'ses-02': 5,
        'ses-03': 5
    },
}

def validate_database():
    """
    Validate the structure of RUNS_DATABASE
    
    Raises:
        ValueError: If database structure is invalid
    """
    if not isinstance(RUNS_DATABASE, dict):
        raise ValueError("RUNS_DATABASE must be a dictionary")
    
    for sub, sessions in RUNS_DATABASE.items():
        if not isinstance(sessions, dict):
            raise ValueError(f"Subject {sub} must have a dictionary of sessions")
        
        for ses, n_runs in sessions.items():
            if not isinstance(n_runs, int):
                raise ValueError(f"Number of runs for {sub} {ses} must be an integer")
            if n_runs <= 0:
                raise ValueError(f"Number of runs for {sub} {ses} must be positive")

# Validate database structure on import
validate_database() 