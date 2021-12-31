
#!/bin/bash
function fail {
  echo $1 >&2
  exit 1
}

function retry {
  local n=1
  local max=10
  local delay=1
  while true; do
    "$@" && break || {
      if [[ $n -lt $max ]]; then
        ((n++))
        echo "Command failed. Attempt $n/$max:"
        sleep $delay;
      else
        fail "The command has failed after $n attempts."
      fi
    }
  done
}

time retry bash run_n.sh -c LI -a TFSI -f EC -n 2 -P LI -N S1 -S O2 -R 2 -V 30
