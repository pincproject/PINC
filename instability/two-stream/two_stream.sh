#!/bin/bash

if [ -d "data" ]; then
  echo "Data folder already exists. Removing old folder..."
  rm -rf data
fi

echo "confirming whether the data folder deleted.." 
ls -l


echo "Creating new data folder..."
mkdir data

echo "Running the program..."
mpirun -np 4 pinc two_stream.ini

echo "Program execution complete."
