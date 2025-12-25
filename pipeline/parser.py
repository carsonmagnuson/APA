import argparse

def parse(inputted: str) -> argparse.Namespace:
  """
    Parses input to determine which pipeline to run and how to run it.

    Args:
        inputted: what is the string of flags that followed the call?

    Returns:
        Should return a list of arguments? I'm not sure how this works yet.

  """
  parser = argparse.ArgumentParser(description="A not so simple command line tool for running an optimized positive selection analysis on a gene of interest.")
  
