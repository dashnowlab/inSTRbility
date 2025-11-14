from filelock import FileLock, Timeout

def is_locked(filepath):
    lock = FileLock(filepath + ".lock")
    try:
        # Try to acquire immediately (non-blocking)
        with lock.acquire(timeout=0):
            return False  # lock was free, we acquired and released it
    except Timeout:
        return True  # someone else holds the lock

# Example
if is_locked("../../genome-references/hg38/GRCh38.M3-complex.bed"):
    print("Compression file is locked (someone is writing)")
else:
    print("No lock, safe to proceed")
