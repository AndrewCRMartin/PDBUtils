#!/usr/local/bin/perl

# The line above tells the shell to run this as a Perl program. You
# need to chage the path to specify wherever you have Perl installed.
# You also need to do chmod +x on this file
# Don't remove that first comment line...



# Use this as a flag to keep track of whether we have an open file
# In Perl all variables begin with a $
$ok = 0;

# Step over each record in the input file
# To read from a file one simply encloses the filehandle (which you
# get by "opening" the file) in <>. Perl automatically opens the
# filename specified on the command line (or standard input from
# a pipe) and knows that if you don't specify a filename that you
# are referring to a file called STDIN which is this input from a
# pipe or the filename on the command line. So this is saying read
# each line from the file one at a time.
while(<>)
{
    # If it's a remark record
    # This line does a pattern match. It looks for the pattern between
    # the slashes in the current input record
    if(/REMARK/)
    {
        # Split out the filename as $cnam
        # The split command takes the current line and splits it at
        # any white space. So the first $junk will get REMARK put into
        # it, $cnam will be the filename and the second $junk will be
        # the score.
        ($junk,$cnam,$junk) = split;

        # If we had a file open, close it
        # We are using $ok as a flag to say whether there is an output
        # file open at the moment. If there is, close it
        close OUT if($ok);

        # Assume we will manage to open a file for writing (so we set
        # the flag to 1 to say that we have done so). Then try to open
        # it. The open() command returns TRUE if it opened the file
        # OK. We want to test if it failed so we put a ! in front of
        # the open() to mean "NOT"
        $ok = 1;
        if(!open(OUT,">$cnam"))
        { 
            # Failed to open it, so print an error message and set the
            # $ok flag to say that we don't have an output file open
            print "Unable to open file for writing: $cnam\n"; 
            $ok = 0;
        }
    }

    # If an output file is open, print the current record to it
    print OUT if($ok);
}


# We could close the last output file here, but this will happen
# automatically when the program exits
