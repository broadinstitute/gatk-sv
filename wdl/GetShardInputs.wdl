version 1.0

workflow GetShardInputs {
  input {
    Array[String] all_items
    Int shard_number
    Int items_per_shard
    Int num_items
  }

  scatter (j in range(items_per_shard)) {
    Int idx = shard_number * items_per_shard + j
    if (idx < num_items) {
      String shard_items_ = all_items[idx]
    }
  }

  output {
    Array[String] shard_items = select_all(shard_items_)
  }
}