#!/usr0/local/bin/python

import sys
import getopt
import re


###################################
# Help String
###################################

def report_error( error_string ):
    print "*******************************"
    print error_string
    print "For help use --help or -h"
    print "*******************************"
    sys.exit()
    return

help_string = \
            """
Usage: python findprobe.py [args]
Arguments:
  --seq <sequence>
  --seqfile <sequence file in FASTA format>
  --tcount <number of T's required in a probe (default = 5)>
  --tsep <distance between T's, (default = 10)>
                 This value can also be a range, i.e.: --tsep 10-12
  --gccontent <n1-n2>
                 This is the range of permissible gc content for a
                 probe.  Default is --gccontent .46-.54)
  --probelength <n1-n2>
                 This is the range of permissible probe lengths.
                 Default is --probelength 46-54
  --drawincontext
                 Report each probe in the context of the initial
                 sequence.  This is generally unreadable in a terminal
                 window.  Save the output to an html file and view in
                 a web browser for best results.  To save to a file:
                 \"python findprobe.py [options] > output.html\".
                 When using this option, each probe is printed with
                 its percent GC content and melting temp.  Each probe
                 also has a link to more information about that probe.

The program finds probes with exactly <tcount> T's, with <tsep>
nucleotides between each T.  A range of values can be given for
--tsep, for example: "--tsep 8-11".  The program reports probes whose
length falls within the range <probelength>, and whose GC content
falls within the range <gccontent>.  The program tries to get each
probe's GC content as close to .50 as possible.  Also, the program
will not allow a T to be the last base in a probe.

The program as written might take a while to find probes.  The
running time grows with the size of the range in <tsep>.  Running time
could be improved in a future version.

Benjamin Vernot
Version: .8
"""

###################################
# Parameters object!
###################################

class parameters:
    t_count = 5
    t_separation_range = (10,)
    gc_content_range = (.46, .54)
    probe_length_range = (46, 54)
    draw_in_context = False


###################################
# Main method
###################################

def main():


    ## parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["seq=",
                                                       "seqfile=",
                                                       "tcount=",
                                                       "gccontent=",
                                                       "tsep=",
                                                       "probelength=",
                                                       "drawincontext",
                                                       "help"])
        pass

    ## throw an error if an incorrect argument is given
    except getopt.error, msg:
        report_error( msg )
        pass

    # process options
    options = {}
    for o, a in opts:
        
        ## print help string
        if o in ("-h", "--help"):
            print help_string
            sys.exit()
        
        options[o] = a

        pass


    ## make sure that either a sequence or a sequence file is given
    if not ( options.has_key( "--seqfile" ) or options.has_key( "--seq" ) ):
        report_error( "The option --seqfile or --seq must be given." )
        pass

    ## the sequence!
    seq = ""

    ## get report style
    if options.has_key( "--drawincontext" ):
        parameters.draw_in_context = True
        pass


    ## get the sequence from the command line
    if options.has_key( "--seq" ):
        seq = clean_seq( options["--seq"] )
        pass


    ## get the sequence from a FASTA file
    if options.has_key( "--seqfile" ):
        seq = read_fasta_file( options["--seqfile"] )
        pass


    ## get parameters.t_count
    if options.has_key( "--tcount" ):
        parameters.t_count = int(options["--tcount"])
        pass

    ## function to extract a range from a variable
    def extract_range( key, type, single_value_ok ):

        ## single value given
        if single_value_ok and len(options[key].split("-")) == 1:
            l = (type(options[key]),)
            pass
        
        ## this is a proper range
        elif len(options[key].split("-")) == 2:
            s1 = options[key].split("-")[0]
            s2 = options[key].split("-")[1]
            l = (type(s1), type(s2))
            pass
        
        ## incorrect number of values given!
        else:
            report_error( "The option %s has an improper argument: %s" % (key, options[key]) )
            pass
        
        return l


    ## get t_separation
    if options.has_key( "--tsep" ):
        l = extract_range( "--tsep", int, single_value_ok = True )
        if len(l) == 1:
            # add one, so that tsep is the bases between T's
            parameters.t_separation_range = l
            pass
        else:
            # add one to both the beginning and end,
            # so that tsep is the bases between T's
            parameters.t_separation_range = range(l[0], l[1]+1)
            pass
        pass

    
    ## get probe length range
    if options.has_key( "--probelength" ):
        parameters.probe_length_range = extract_range( "--probelength", int, single_value_ok = False )
        pass


    ## get gc content range
    if options.has_key( "--gccontent" ):
        parameters.gc_content_range = extract_range( "--gccontent", float, single_value_ok = False )
        pass


    ## report standard info
    if parameters.draw_in_context: print "<pre>"
    print "Sequence length: %d" % len( seq )
    print "Number of T's required: %d" % parameters.t_count
    print "Distance between each T: %s" % parameters.t_separation_range
    print "Probe length range: %d-%d" % parameters.probe_length_range
    print "GC content range: %.2f-%.2f" % parameters.gc_content_range

    ## initialize the gc content arrays
    # print "... initializing gc count arrays ... "
    init_gc_count_table( seq )

    ## find probes
    # print "... finding probes ... "
    probes = find_variable_tsep_probes( seq = seq )

    ## sort probes based on starting point
    def sort_fn( a, b ):
        if a.start != b.start:
            return a.start - b.start
        if a.t_spots[0] != b.t_spots[0]:
            return a.t_spots[0] - b.t_spots[0]
        if a.t_spots[1] != b.t_spots[1]:
            return a.t_spots[1] - b.t_spots[1]
        return 0
    probes.sort( sort_fn )
    
    ## report probes
    report_probes( seq, probes )
    
    return




###################################
# Read FASTA file
###################################

def read_fasta_file( fasta_file ):

    seq = ""
    
    f = open( fasta_file, "U" )
    try:
        for line in f:
            if len(line) > 0 and line[0] != '>':
                seq += line
                pass
            pass
        pass
    finally:
        f.close()
        pass
    
    return clean_seq( seq )


###################################
# Clean the sequence
###################################

def clean_seq( seq ):
    ## remove any whitespace
    return re.sub("\s+", "", seq)


###################################
# Convert a base to its RC
###################################

def convert(c):
    if c == 'T':
        return 'A'
    if c == 'A':
        return 'T'
    if c == 'G':
        return 'C'
    if c == 'C':
        return 'G'
    pass


###################################
# Initialize the gc count table in parameters
###################################

def init_gc_count_table( seq ):
    
    # the outside array
    parameters.gc_count_table = []

    # list of potential probe sizes
    probe_size_list = range( parameters.probe_length_range[0],
                              parameters.probe_length_range[1] + 1 )
    
    # the inside arrays, one for each window
    for window_size in probe_size_list:
        
        # in these arrays, one spot for each potential starting spot
        parameters.gc_count_table.append( range( len(seq) - window_size + 1 ) )
        
        pass
    
    # fill the things
    for (i, window_array) in enumerate( parameters.gc_count_table ):
        
        window_size = probe_size_list[i]
        
        # set up the initial probe
        p = probe()
        p.start = 0
        p.end = window_size
        p.calculate_gc_content( seq )
        
        for pos in xrange( len(seq) - window_size + 1 ):
            
            p.calculate_gc_content_for_range( seq, pos, pos + window_size )
            parameters.gc_count_table[i][pos] = p.gc_count
            
            pass
        
        pass
    
    return



###################################
# Find the probes allowing variable tseps
###################################

def find_variable_tsep_probes( seq ):

    ## debug variable
    debug = False
    
    ## list of completed probes
    completed_probes = []

    ## hash table
    ## the keys are sequence positions where we are interested in seeing a T
    ## the values are a list of probes which would "like" to have a T at this position.
    ## - these probes are not yet completed, and are mostly used for their t_spots list
    potential_probes = {}

    ## function which modifies the potential_probes list
    def populate_potential_probes( t_spots, index, print_string = None ):

        if print_string and debug: print print_string
        
        ## add a probe to the list of potential probes,
        ## once for each value in tseprange
        for tsep in parameters.t_separation_range:
            
            ## make new probes based on the original probe
            p = probe()
            p.t_spots = list(t_spots)
            p.start = t_spots[0]
            p.end = t_spots[-1] + 1

            ## the next index we would look for a T
            new_index = index + tsep + 1
            if print_string and debug:
                print "%d: %s: %s: %d" % (p.start, p.tspots_string( seq ), p.t_spots, new_index)
                pass
            
            if potential_probes.has_key( new_index ):
                potential_probes[new_index].append( p )
                pass
            else:
                potential_probes[new_index] = [ p ]
                pass
            pass
        return
    

    for index, nucleotide in enumerate(seq):

        if debug: print "at nucleotide %s, index %d" % (nucleotide, index)

        ## this is a T!
        if nucleotide == 'A':

            ####
            ## create a potential probe which starts at this T
            populate_potential_probes( t_spots = [index],
                                        index = index,
                                        print_string = "potential probe" )


            ####
            ## now check for probes we are already building

            ## we were expecting a T at this location!
            if potential_probes.has_key( index ):


                ## get a list of each probe which expected a T in this position
                probe_list = potential_probes[index]


                ## for each probe, check to see if this T completes it
                for p in list(probe_list):

                    ## this T completes this probe
                    if len(p.t_spots) + 1 == parameters.t_count:
                        
                        ## add this T
                        p.t_spots += [index]
                        
                        ## update the stats for this probe
                        p.compute_stats( seq )
                        
                        ## if the gc_content is correct, add it to the list of completed probes
                        if p.gc_content:
                            completed_probes.append( p )
                            if debug: print "found a good probe: " + p.tspots_string( seq )
                            pass
                        else:
                            if debug: print "removing a probe: " + p.tspots_string( seq )
                            pass

                        ## remove it from the list of potential probes, even if we did not
                        ## add it to the list of completed probes
                        probe_list.remove( p )
   
                        pass

                    else:
                        if debug: print "not considering probe: " + p.tspots_string( seq )
                    
                    pass
                
                ## for each probe, we have to mark the next spot we would expect to see a T
                ## we also keep track of where this probe has seen T's
                for p in probe_list:
                    populate_potential_probes( t_spots = p.t_spots + [index],
                                                index = index,
                                                print_string = "growing probe at %d" % index )
                    pass

                pass

            pass
        
        ## remove the probes which expected T's here.
        ## so sad!
        if potential_probes.has_key( index ):
            del potential_probes[ index ]
            pass


#        if debug: print "  after  %s" % potential_probes

        # end of loop through sequence #
        pass
    
    
    return completed_probes





###################################
# Probe Report
###################################

def report_probes( seq, probes ):

    last_bp = 0

    print "%d probes found" % len( probes )
    
    for (p_i, p) in enumerate( probes ):

        if last_bp <= p.start and last_bp != 0:
            print ""
            print "|-------------------------|"
            print "|-- gap between probes  --| %d bp" % (p.start - last_bp)
            print "|-------------------------|"
            pass
        
        print ""
        if parameters.draw_in_context:
            print "<a name=\"%d\" href=\"#%ddic\">Probe #%d</a>" % (p_i, p_i, p_i)
            pass
        else:
            print "Probe #%d" % p_i
            pass
        print "Start of probe: %d" % p.start
        print "End of probe: %d" % p.end
        print "Number of base pairs: %d" % (p.end - p.start)
        print "Percentage GC content: %.2f" % p.gc_content
        print "Melting point: %2.2f*C" % p.melting_point

        ## print the probe
        print "Sequence:"
        print seq[p.start:p.end]
        ## print the tspots
        print p.tspots_string( seq, "forward" )
        
        ## print the reverse compliment
        print "Probe (reverse complement):"
        print p.reverse_compliment( seq )
        ## print the tspots
        print p.tspots_string( seq, "reverse" )
        ## print in the alternate format
        print "Probe (reverse complement, alternate format):"
        print p.reverse_compliment_probe_formatted( seq )

        ## keep track of the last place on the sequence that we've seen
        if last_bp < p.end:
            last_bp = p.end
            pass
        
        pass

    print ""
    print "%d probes found" % len( probes )

    if parameters.draw_in_context:

        leading_space = len("Query Sequence: ")
        
        print ""
        print ""
        print "Query Sequence: " + seq

        for (p_i, p) in enumerate( probes ):

            pstring = ""
            
            # print leading space
            for i in xrange(0, p.start + leading_space -
                            len("Probe #%d: " % p_i) ):
                pstring += " "
                pass
            pstring += "<a href=\"#%d\">Probe #%d:</a> <a name=\"%ddic\">" % (p_i, p_i, p_i)

            # print sequence
            pstring += seq[p.start:p.end]

            # print temp range
            pstring += " | %2.2f degrees Celsius     " % p.melting_point

            # close the anchor tag
            pstring += "</a>"
            
            # actually print
            print pstring

            pstring = ""
            
            # print leading space
            for i in xrange(0, p.start + leading_space):
                pstring += " "
                pass

            # print tspots
            pstring += p.tspots_string( seq )

            # print gc content range
            pstring += " | %.2f %% GC content" % p.gc_content

            # actually print
            print pstring

    return




###########################################
# probe class
###########################################

class probe:

    def __init__(self):
        self.gc_content = None
        self.gc_count = None
        self.start = None
        self.end = None
        self.t_spots = None
        self.melting_point = None

        self.gc_min = 1
        self.gc_max = 0
        self.gc_deviation = 1
        self.temp_min = sys.maxint
        self.temp_max = 0
        return

    def __repr__(self):
        return "p:%s" % self.t_spots

    ## update the stats for this probe
    ## based entirely on contents of t_spots
    def compute_stats( self, seq ):

        debug = False

        ## keep track of min and max temp and gc
        self.gc_min = 1
        self.gc_max = 0
        self.gc_deviation = 1
        self.temp_min = sys.maxint
        self.temp_max = 0

        ## save a choice probe
        best_start = None
        best_end = None

        ## make sure that the start and end are correct
        self.start = self.t_spots[0]
        self.end = self.t_spots[-1] + 1

        ## get the length of the probe from T to T
        tlength = self.end - self.start
        
        ## window sizes we have to consider
        window_sizes = range( max(tlength, parameters.probe_length_range[0]),
                              max(tlength, parameters.probe_length_range[1] + 1) )

        ## consider each window size
        for window_size in window_sizes:

            ## get bounds on the starting points we have to consider
            s = max( 0, self.t_spots[0] - window_size + tlength )
            # don't allow the window to go all the way up against the first A
            e = min( len(seq), self.t_spots[0] + window_size ) - window_size
            if debug:
                print "window size: %d" % window_size
                print "tlength: %d" % tlength
                print "first tspot: %d" % self.t_spots[0]
                pass
            
            ## loop through each placement of the window
            for starting_point in xrange(s, e):

                ## get the gc content for this window
                new_start = starting_point
                new_end = starting_point + window_size
                # this function also sets the new start and end for self
                self.fetch_gc_content_for_range( seq, new_start, new_end )

                ## print debug data
                if debug:
                    ## print the original probe from t to t
                    debug_string = "".join( [" " for i in xrange(0, self.t_spots[0])] )
                    debug_string += seq[self.t_spots[0]:self.t_spots[-1]+1]
                    print debug_string
                    ## print the tspots
                    debug_string = "".join( [" " for i in xrange(0, self.start)] )
                    debug_string += self.tspots_string( seq )
                    print debug_string
                    ## print the window we're currently considering
                    debug_string = "".join( [" " for i in xrange(0, self.start)] )
                    debug_string += "|"
                    debug_string += "".join( ["." for i in xrange(self.start+1, self.end-1)] )
                    debug_string += "| " + str(self.gc_content)
                    print debug_string
                    pass

                ## calculate the melting point
                self.melting_point = 64.9 + 41 * (self.gc_count - 16.4) / (self.end - self.start)
                
                ## see if this probe passes the gc content test
                if self.gc_content >= parameters.gc_content_range[0] and \
                       self.gc_content <= parameters.gc_content_range[1] :

                    if abs( self.gc_content - .5 ) < self.gc_deviation:
                        self.gc_deviation = abs( self.gc_content - .5 )
                        best_start = self.start
                        best_end = self.end
                        pass
                    
                    if self.gc_content < self.gc_min:
                        self.gc_min = self.gc_content
                        pass
                    
                    if self.gc_content > self.gc_max:
                        self.gc_max = self.gc_content
                        pass
                    
                    if self.melting_point < self.temp_min:
                        self.temp_min = self.melting_point
                        pass
                    
                    if self.melting_point > self.temp_max:
                        self.temp_max = self.melting_point
                        pass

                    # end if statement checking gc content
                    pass
                # end for loop moving this window from left to right
                pass

            ## take the shortest sequence with the highest gc content
            if False and best_start:
                break
            
            # end for loop checking each window
            pass

        ## reset the values for this probe and return
        if best_start:
            self.start = best_start
            self.end = best_end
            self.calculate_gc_content( seq )
            self.melting_point = 64.9 + 41 * (self.gc_count - 16.4) / (self.end - self.start)
            pass
        else:
            self.gc_content = None
            self.gc_count = None
            self.melting_point = None
            self.start = self.t_spots[0]
            self.end = self.t_spots[-1] + 1
            pass

        return

    ## function to print tspots, including "*" under each important T
    ## and the number of spots between them
    def tspots_string( self, seq, orientation = "forward" ):
        tspots = ""
        for i, t in enumerate( self.t_spots ):
            if i == 0:
                tspots = "".join( [' ' for ii in range(self.start, t)] )
                pass
            if i != len( self.t_spots ) - 1:
                dist = "*   "
                if orientation == "forward": dist += str( self.t_spots[i+1] - self.t_spots[i] - 1 )
                if orientation == "reverse": dist += str( self.t_spots[i+1] - self.t_spots[i] - 1 )[::-1]
                t2 = self.t_spots[i+1] - len(dist)
                tspots += dist
                tspots += "".join( [' ' for ii in range(t, t2)] )
                pass
            else:
                tspots += "*" + "".join( [' ' for ii in range(t, self.end - 1)] )
                pass
            
            pass

        if orientation == "reverse":
            tspots = tspots[::-1]
            pass
        elif orientation != "forward":
            report_error( "Bad orientation in tspots_string: %s" % orientation )
            pass
        
        return tspots
    
    ## functions to get teh reverse compliment of this probe
    def reverse_compliment( self, seq ):
        return "".join ([convert(c)
                         for c in seq[self.end-1:self.start-1:-1]])
    
    def reverse_compliment_probe_formatted( self, seq ):
        rc = ""

        for e, i in enumerate( xrange(self.end-1, self.start-1, -1) ):

            if i in self.t_spots:
                rc += "*%s*" % convert(seq[i])
                pass
            else:
                rc += convert(seq[i])
                pass

            if ((e+1) % 3) == 0:
                rc += " "
                pass
            pass

        return rc
            
            
    
    ## set the percentage gc content and gc count for a probe
    def calculate_gc_content( self, seq ):
        
        self.gc_count = 0
        
        ## look at every nucleotide from start to end of the sequence
        for index, nucleotide in enumerate( seq[self.start:self.end] ):
            
            ## is a G or C
            if nucleotide == 'G' or nucleotide == 'C':
                self.gc_count += 1
                pass
            
            pass
        
        self.gc_content = float( self.gc_count ) / float( self.end - self.start )
        
        return
    

    ## set the percentage gc content and gc count for a probe
    def calculate_gc_content_for_range( self, seq, start, end ):
        
        ## look at every nucleotide from new start to old start and add new G or C
        for index, nucleotide in enumerate( seq[start:self.start] ):
            if nucleotide == 'G' or nucleotide == 'C': self.gc_count += 1
            pass
        
        ## look at every nucleotide from old start to new start and subtract new G or C
        for index, nucleotide in enumerate( seq[self.start:start] ):
            if nucleotide == 'G' or nucleotide == 'C': self.gc_count -= 1
            pass
        
        ## look at every nucleotide from old end to new end and add new G or C
        for index, nucleotide in enumerate( seq[self.end:end] ):
            if nucleotide == 'G' or nucleotide == 'C': self.gc_count += 1
            pass
        
        ## look at every nucleotide from old end to new end and subtract new G or C
        for index, nucleotide in enumerate( seq[end:self.end] ):
            if nucleotide == 'G' or nucleotide == 'C': self.gc_count -= 1
            pass
        
        self.gc_content = float( self.gc_count ) / float( end - start )

        self.start = start
        self.end = end
        
        return

    
    ## set the percentage gc content and gc count for a probe
    def fetch_gc_content_for_range( self, seq, start, end ):

        window_pos = end - start - parameters.probe_length_range[0]
        
        self.gc_count = parameters.gc_count_table[window_pos][start]
        self.start = start
        self.end = end
        self.gc_content = float( self.gc_count ) / float( self.end - self.start )
        
        return
    
#### end of probe class
    pass

    


#########################################################################################################
#########################################################################################################
# Invoke the main method
#########################################################################################################
#########################################################################################################

if __name__ == "__main__":
    main()
    pass

#########################################################################################################
#########################################################################################################
# Invoke the main method
#########################################################################################################
#########################################################################################################

