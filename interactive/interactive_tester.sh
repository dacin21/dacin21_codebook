#!/bin/bash

TMPDIR="/home/dacin21/tmp"

if [[ "$1" = "-c" ]]; then
  echo "Cleaning up tmp dir"
  rm -r "${TMPDIR}"/*
  exit 0
fi

FORCE=0
if [[ "$1" = "-f" ]]; then
  FORCE=1
  shift
fi

if [[ -z "$2" ]]; then
  echo "usage: $0 [-c] [-f] <judge> <solution> [compile command]"
  exit 1
fi

JUDGE="$1"
SOL="$2"

if ! [[ -z "$3" ]]; then
echo "compiling solution"
 if ! $($3); then
   exit 1
 fi
fi

# interupt all panes
tmux send-keys -t sol_stderr C-c
tmux send-keys -t judge_stderr C-c
tmux send-keys -t judge C-c
tmux send-keys -t sol C-c

# setup fifo files
MYTMPDIR=$(mktemp -d --tmpdir=${TMPDIR})
if [[ -z ${MYTMPDIR} ]]; then
  echo "failed to create tmp dir"
  exit 1
fi
# trap "rm -rf $MYTMPDIR" EXIT # TODO: find a way to clean up tmp dir
mkfifo -m 600 "$MYTMPDIR/sol_out"
mkfifo -m 600 "$MYTMPDIR/sol_err"
mkfifo -m 600 "$MYTMPDIR/judge_out"
mkfifo -m 600 "$MYTMPDIR/judge_err"

tmux send-keys -t sol_stderr "cat $MYTMPDIR/sol_err" C-m
tmux send-keys -t judge_stderr "cat $MYTMPDIR/judge_err" C-m
tmux send-keys -t judge "${JUDGE} < \"$MYTMPDIR/sol_out\" 2> \"$MYTMPDIR/judge_err\" | tee \"$MYTMPDIR/judge_out\"" C-m
tmux send-keys -t sol "${SOL} < \"$MYTMPDIR/judge_out\" 2> \"$MYTMPDIR/sol_err\" | tee \"$MYTMPDIR/sol_out\"" C-m
