#!/usr/bin/env python3
"""
CLI script to copy Google Cloud Storage objects from a JSON file.
Processes JSON where keys can be GCS URIs (or lists of URIs) and copies files
from buckets matching a specified prefix to a destination bucket.
"""

import argparse
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, Optional, Tuple, Union, List


def parse_gcs_uri(uri: str) -> tuple[str, str]:
    """
    Parse a GCS URI into bucket and path components.
    
    Args:
        uri: GCS URI in format gs://bucket/path/to/file
        
    Returns:
        Tuple of (bucket_name, object_path)
    """
    match = re.match(r'gs://([^/]+)/(.+)', uri)
    if not match:
        raise ValueError(f"Invalid GCS URI: {uri}")
    return match.group(1), match.group(2)


def should_copy_bucket(bucket_name: str, source_buckets: List[str]) -> bool:
    """
    Check if bucket name matches one of the specified source buckets.
    
    Args:
        bucket_name: Name of the GCS bucket
        source_buckets: List of bucket names to match exactly
        
    Returns:
        True if bucket should be copied
    """
    return bucket_name in source_buckets


def resolve_tasks_in_value(value: Any, task_results: Dict[int, Optional[str]]) -> Any:
    """
    Resolve task placeholders in a value to actual URIs based on task results.
    
    Args:
        value: Value with task placeholders
        task_results: Dict mapping task index to destination URI (or None on failure)
        
    Returns:
        Value with task placeholders replaced by actual URIs
    """
    if isinstance(value, dict):
        if value.get("type") == "literal":
            return value["value"]
        elif value.get("type") == "task":
            task_index = value["index"]
            result_uri = task_results.get(task_index)
            # If copy failed, use original URI
            return result_uri if result_uri else value["original"]
        else:
            return value
    elif isinstance(value, list):
        return [resolve_tasks_in_value(item, task_results) for item in value]
    else:
        return value


def execute_copy_tasks(copy_tasks: List[Dict[str, Any]], dry_run: bool, max_workers: int) -> Tuple[Dict[int, Optional[str]], Dict[str, int]]:
    """
    Execute copy tasks in parallel using ThreadPoolExecutor.
    
    Args:
        copy_tasks: List of copy task dictionaries
        dry_run: If True, simulate copies without executing
        max_workers: Maximum number of parallel workers
        
    Returns:
        Tuple of (task_results, stats)
        - task_results: Dict mapping task index to destination URI (or None on failure)
        - stats: Statistics about copy operations
    """
    task_results = {}
    stats = {"copied": 0, "errors": 0}
    
    if not copy_tasks:
        return task_results, stats
    
    print(f"\nExecuting {len(copy_tasks)} copy tasks with {max_workers} parallel workers...")
    print("-" * 60)
    
    def copy_task(task):
        source_uri = task["source_uri"]
        dest_bucket = task["dest_bucket"]
        task_index = task["task_index"]
        
        dest_uri = copy_gcs_file(source_uri, dest_bucket, dry_run)
        return task_index, dest_uri
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_task = {executor.submit(copy_task, task): task for task in copy_tasks}
        
        # Process completed tasks
        completed = 0
        for future in as_completed(future_to_task):
            completed += 1
            try:
                task_index, dest_uri = future.result()
                task_results[task_index] = dest_uri
                
                if dest_uri:
                    stats["copied"] += 1
                else:
                    stats["errors"] += 1
                    
                # Print progress
                if completed % 10 == 0 or completed == len(copy_tasks):
                    print(f"Progress: {completed}/{len(copy_tasks)} files processed")
                    
            except Exception as e:
                task = future_to_task[future]
                print(f"Error executing task for {task['source_uri']}: {e}", file=sys.stderr)
                task_results[task["task_index"]] = None
                stats["errors"] += 1
    
    print("-" * 60)
    return task_results, stats


def copy_gcs_file(source_uri: str, dest_bucket: str, dry_run: bool = False) -> Optional[str]:
    """
    Copy a GCS file to the destination bucket, maintaining subdirectory structure.
    
    Args:
        source_uri: Source GCS URI
        dest_bucket: Destination bucket name (without gs://)
        dry_run: If True, only print what would be copied
        
    Returns:
        Destination URI if successful, None otherwise
    """
    try:
        source_bucket, object_path = parse_gcs_uri(source_uri)
        dest_uri = f"gs://{dest_bucket}/{object_path}"
        
        if dry_run:
            print(f"[DRY RUN] Would copy: {source_uri} -> {dest_uri}")
            return dest_uri
        
        print(f"Copying: {source_uri} -> {dest_uri}")
        result = subprocess.run(
            ["gsutil", "cp", source_uri, dest_uri],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"Error copying {source_uri}: {result.stderr}", file=sys.stderr)
            return None
        
        return dest_uri
        
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Unexpected error copying {source_uri}: {e}", file=sys.stderr)
        return None


def process_json_value(value: Any, source_buckets: List[str], dest_bucket: str, copy_tasks: List[Dict[str, Any]]) -> Tuple[Any, Dict[str, int]]:
    """
    Process a JSON value that could be a URI string or list of URIs.
    Collects copy tasks instead of executing them immediately.
    
    Args:
        value: JSON value (string or list)
        source_buckets: List of source bucket names to filter on
        dest_bucket: Destination bucket
        copy_tasks: List to append copy tasks to
        
    Returns:
        Tuple of (updated_value, stats_dict)
        - updated_value: Placeholder with task indices for later replacement
        - stats_dict: Dictionary with counts of processed, copied, and skipped files
    """
    stats = {"processed": 0, "copied": 0, "skipped": 0, "errors": 0}
    
    # Handle both single URI and list of URIs
    if isinstance(value, str):
        is_single = True
        uris = [value]
    elif isinstance(value, list):
        is_single = False
        uris = value
    else:
        # Not a URI or list of URIs, return as-is
        return value, stats
    
    result_uris = []
    
    for uri in uris:
        if not isinstance(uri, str):
            result_uris.append({"type": "literal", "value": uri})
            continue
            
        # Check if it's a GCS URI
        if not uri.startswith("gs://"):
            result_uris.append({"type": "literal", "value": uri})
            continue
        
        stats["processed"] += 1
        
        try:
            bucket_name, _ = parse_gcs_uri(uri)
            
            if should_copy_bucket(bucket_name, source_buckets):
                # Add to copy tasks list
                task_index = len(copy_tasks)
                copy_tasks.append({
                    "source_uri": uri,
                    "dest_bucket": dest_bucket,
                    "task_index": task_index
                })
                result_uris.append({"type": "task", "index": task_index, "original": uri})
            else:
                stats["skipped"] += 1
                result_uris.append({"type": "literal", "value": uri})
                
        except Exception as e:
            print(f"Error processing URI {uri}: {e}", file=sys.stderr)
            stats["errors"] += 1
            result_uris.append({"type": "literal", "value": uri})
    
    # Return single item or list based on input type
    result = result_uris[0] if is_single else result_uris
    return result, stats


def process_json_file(json_file: str, source_buckets: List[str], dest_bucket: str, output_file: Optional[str] = None, dry_run: bool = False, max_workers: int = 8) -> None:
    """
    Process a JSON file and copy matching GCS objects.
    
    Args:
        json_file: Path to JSON file
        source_buckets: List of source bucket names to filter on
        dest_bucket: Destination bucket name
        output_file: Path to output JSON file with updated URIs
        dry_run: If True, only print what would be copied
        max_workers: Maximum number of parallel copy workers
    """
    print(f"Processing JSON file: {json_file}")
    print(f"Source buckets: {', '.join(source_buckets)}")
    print(f"Destination bucket: {dest_bucket}")
    if output_file:
        print(f"Output JSON file: {output_file}")
    if dry_run:
        print("DRY RUN MODE - No files will be copied")
    print("-" * 60)
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File not found: {json_file}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {json_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # First pass: collect all copy tasks
    copy_tasks = []
    total_stats = {"processed": 0, "copied": 0, "skipped": 0, "errors": 0}
    placeholder_data = None
    
    # Process the JSON structure to collect tasks
    if isinstance(data, dict):
        placeholder_data = {}
        for key, value in data.items():
            result_value, stats = process_json_value(value, source_buckets, dest_bucket, copy_tasks)
            placeholder_data[key] = result_value
            total_stats["processed"] += stats["processed"]
            total_stats["skipped"] += stats["skipped"]
            total_stats["errors"] += stats["errors"]
    elif isinstance(data, list):
        placeholder_data = []
        for value in data:
            result_value, stats = process_json_value(value, source_buckets, dest_bucket, copy_tasks)
            placeholder_data.append(result_value)
            total_stats["processed"] += stats["processed"]
            total_stats["skipped"] += stats["skipped"]
            total_stats["errors"] += stats["errors"]
    else:
        placeholder_data = data
    
    # Second pass: execute copy tasks in parallel
    task_results, copy_stats = execute_copy_tasks(copy_tasks, dry_run, max_workers)
    total_stats["copied"] = copy_stats["copied"]
    total_stats["errors"] += copy_stats["errors"]
    
    # Third pass: resolve task placeholders to actual URIs
    updated_data = None
    if isinstance(placeholder_data, dict):
        updated_data = {k: resolve_tasks_in_value(v, task_results) for k, v in placeholder_data.items()}
    elif isinstance(placeholder_data, list):
        updated_data = [resolve_tasks_in_value(v, task_results) for v in placeholder_data]
    else:
        updated_data = placeholder_data
    
    # Write output JSON file if specified
    if output_file and updated_data is not None:
        try:
            with open(output_file, 'w') as f:
                json.dump(updated_data, f, indent=2)
            print("-" * 60)
            print(f"Updated JSON written to: {output_file}")
        except Exception as e:
            print(f"Error writing output file {output_file}: {e}", file=sys.stderr)
    
    # Print summary
    print("-" * 60)
    print(f"Summary:")
    print(f"  URIs processed: {total_stats['processed']}")
    print(f"  Files copied: {total_stats['copied']}")
    print(f"  Files skipped (bucket not in source list): {total_stats['skipped']}")
    print(f"  Errors: {total_stats['errors']}")


def main():
    parser = argparse.ArgumentParser(
        description="Copy GCS objects from JSON file to destination bucket",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s input.json my-dest-bucket --source-buckets fc-bucket1,fc-bucket2

The JSON file should have GCS URIs as keys (or lists of URIs).
Only files from the specified source buckets will be copied.
        """
    )
    
    parser.add_argument(
        "json_file",
        help="Path to JSON file containing GCS URIs"
    )
    
    parser.add_argument(
        "dest_bucket",
        help="Destination GCS bucket name (without gs://)"
    )
    
    parser.add_argument(
        "--source-buckets",
        "-s",
        required=True,
        help="Comma-separated list of source bucket names to copy from"
    )
    
    parser.add_argument(
        "--output",
        "-o",
        help="Path to output JSON file with updated URIs (default: <input>_updated.json)"
    )
    
    parser.add_argument(
        "--workers",
        "-w",
        type=int,
        default=8,
        help="Number of parallel workers for copying files (default: 8)"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be copied without actually copying"
    )
    
    args = parser.parse_args()
    
    # Parse comma-separated source buckets
    source_buckets = [b.strip() for b in args.source_buckets.split(',')]
    # Remove gs:// prefix if present
    source_buckets = [b[5:] if b.startswith("gs://") else b for b in source_buckets]
    
    # Validate destination bucket format
    if args.dest_bucket.startswith("gs://"):
        args.dest_bucket = args.dest_bucket[5:]
    
    # Determine output file
    output_file = args.output
    if output_file is None:
        # Generate default output filename
        import os
        base, ext = os.path.splitext(args.json_file)
        output_file = f"{base}_updated{ext}"
    
    process_json_file(args.json_file, source_buckets, args.dest_bucket, output_file, args.dry_run, args.workers)


if __name__ == "__main__":
    main()
