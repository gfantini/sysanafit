#!/usr/bin/perl

#(v03) [UNIFIED VERSION] write configfile given 1st argument (allowed both numbers 1..N and strings!)

open($configfile,'>',$ARGV[0]."/".$ARGV[0].".config");

@nsigma = (-2.0,-1.0,1.0,2.0);
for($i=0;$i<10;$i++)
{
    if( $ARGV[0] eq "default") # eq operator confronts STRINGWISE
    {
        print $configfile "0\n";
    }else{
        $resto = $ARGV[0]%4;
        $var = ($ARGV[0]-$resto)/4;
        
        if($i == $var){
            print $configfile "$nsigma[$resto]\n";
        }else{
            print $configfile "0\n";
        }
    }
}
close($configfile);
