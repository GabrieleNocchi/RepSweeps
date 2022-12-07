#!/usr/bin/perl


# Open the VCF file
my $filename = $ARGV[0];           # store the 1st argument into the variable
open $fh, '<', $filename or die $!;

# Initialize variables
$current_chrom = "";
$first_pos = 0;
$last_pos = 0;

# Loop through each line of the VCF file
while(<$fh>) {
  # Ignore comment lines
  if ($_ =~ /^#/) {
    next;
  }

  # Parse the line
  chomp;
  ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $samples) = split("\t", $_);

  # If this is a new chromosome, print the first and last position of the previous chromosome
  if ($chrom ne $current_chrom && $current_chrom ne "") {
    print "$current_chrom\t$first_pos\t$last_pos\n";
    $first_pos = 0;
    $last_pos = 0;
  }

  # Update the current chromosome and the first and last positions
  $current_chrom = $chrom;
  $first_pos = $first_pos == 0 ? $pos : $first_pos;
  $last_pos = $pos;
}

# Print the first and last position of the last chromosome
print "$current_chrom\t$first_pos\t$last_pos\n";

# Close the VCF file
close(VCF);

