#!/bin/bash

FORCE=0
if [[ "$1" = "-f" ]]; then
  FORCE=1
  shift
fi

if [[ -z "$1" ]]; then
  echo "usage: $0 [-f] <session_name>"
  exit 1
fi

SESSION="$1"

# setup tmux session
if tmux has-session -t $SESSION; then
  echo "Session '$SESSION' already exists"
  if [[ $FORCE -eq 1 ]]; then
    echo "Kill session"
    tmux kill-session -t $SESSION
  fi
fi
if ! $(tmux has-session -t $SESSION); then
  echo "Create session"
  tmux new-session -d -s $SESSION
  tmux new-window -t $SESSION:1 -n 'judge'
  tmux new-window -t $SESSION:2 -n 'sol'
  tmux new-window -t $SESSION:3 -n 'judge_stderr'
  tmux new-window -t $SESSION:4 -n 'sol_stderr'
fi

tmux attach-session -t $SESSION
