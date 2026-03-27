import argparse
import re

parser = argparse.ArgumentParser(
    "Use regex to match a pattern in the read ID format this into our expected read IDs"
)
parser.add_argument(
    "regex",
    help="regex of pattern to match on read ID lines. output will be all grouped items. all lines must match",
)
parser.add_argument("input", help="fastq file to read from")
parser.add_argument("output", help="out file to write to")

args = parser.parse_args()

regex = re.compile(args.regex)

with open(args.input, "rt") as input_fastq:
    with open(args.output, "wt") as output_fastq:
        i = 1
        while True:
            id_line = input_fastq.readline()
            data_line = input_fastq.readline()
            separator_line = input_fastq.readline()
            quality_line = input_fastq.readline()

            if id_line == "":
                break  # no more to read

            m = regex.match(id_line)
            if m:
                new_id_line = "".join(m.groups())
            else:
                raise Exception(
                    f"Regex {regex} did not match ID line on line {i}:\n{id_line}"
                )

            output_fastq.write(new_id_line)
            output_fastq.write("\n")
            output_fastq.write(data_line)
            output_fastq.write(separator_line)
            output_fastq.write(quality_line)

            i += 4
