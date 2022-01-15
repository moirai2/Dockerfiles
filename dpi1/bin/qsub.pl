#!/usr/bin/env perl
use strict 'vars';
use Cwd;
use File::Basename;
use File::Temp qw/tempfile tempdir/;
use FileHandle;
use Time::Local;
use Time::localtime;
############################## HEADER ##############################
my($program_name,$prgdir,$program_suffix)=fileparse($0);
$prgdir=substr($prgdir,0,-1);
my $program_version="2021/12/23";
############################## MAIN ##############################
if(scalar(@ARGV)==1&&$ARGV[0]eq"sort"){sortSubs();exit();}
my $rootdir=`pwd`;
chomp($rootdir);
my $qsubdir="$rootdir/qsub";
mkdir($qsubdir);
my $bashdir="$qsubdir/bash";
mkdir($bashdir);
my $stdoutdir="$qsubdir/stdout";
mkdir($stdoutdir);
my $stderrdir="$qsubdir/stderr";
mkdir($stderrdir);
my $quietMode;
if(scalar(@ARGV)==1&&$ARGV[0]eq"test"){test();exit();}
my $options={};
my $arguments=[];
for(my $i=0;$i<scalar(@ARGV);$i++){
    my $arg=$ARGV[$i];
    if($arg eq "-cwd"){$options->{$arg}="use";}
    elsif($arg eq "-v"){
        if(!defined($options->{$arg})){$options->{$arg}={};}
        my $string=$ARGV[++$i];
        my @tokens=split(/\,/,$string);
        foreach my $token(@tokens){
            my ($key,$val)=split(/\=/,$token);
            $options->{$arg}->{$key}=$val;
        }
    }elsif($arg eq "-o"){$options->{$arg}=$ARGV[++$i];}
    elsif($arg eq "-e"){$options->{$arg}=$ARGV[++$i];}
    elsif($arg eq "-N"){$options->{$arg}=$ARGV[++$i];}
    elsif($arg eq "-hold_jid"){$options->{$arg}=$ARGV[++$i];}
    elsif($arg eq "-sync"){$options->{$arg}=$ARGV[++$i];}
    elsif($arg eq "--sleep"){$options->{$arg}=$ARGV[++$i];}
    elsif($arg eq "--maxjob"){$options->{$arg}=$ARGV[++$i];}
    else{push(@{$arguments},$arg)}
}
my $sleepTime=exists($options->{"--sleep"})?$options->{"--sleep"}:60;
my $maxJob=exists($options->{"--maxjob"})?$options->{"--maxjob"}:92;
if(scalar(@{$arguments})==0){exit();}
waitForFreeJobQueue();
my $bashFile=createBash($options,join(" ",@{$arguments}));
my $command="bash $bashFile";
if(!exists($options->{"-sync"})){$command.=" &";}
print STDERR ">$command\n";
system($command);
############################## createBash ##############################
sub createBash{
    my $options=shift();
    my $command=shift();
    if(exists($options->{"-hold_jid"})){
        my $id=$options->{"-hold_jid"};
        my $holdDir="$bashdir/$id";
        my @files=listFilesRecursively("\.sh\$",undef,0,$holdDir);
        if(!quietMode){print STDERR ">Waiting for job '$id' to complete..";}
        while(-e $holdDir&&scalar(@files)>0){
            if(!quietMode){print STDERR ".";}
            @files=listFilesRecursively("\.sh\$",undef,0,$holdDir);
        }
        if(!quietMode){print STDERR "OK\n";}
    }
    my @commands=("#!/bin/bash");
    my ($writer,$tmpFile)=tempfile("bashXXXXXXXXXX",DIR=>$qsubdir,SUFFIX=>".sh");
    my $basename=basename($tmpFile,".sh");
    my $jobDir="$bashdir".(exists($options->{"-N"})?"/".$options->{"-N"}:"");
    if(!-e $jobDir){mkdir($jobDir);}
    my $bashFile="$jobDir/$basename.sh";
    #if(exists($options->{"-cwd"})){push(@commands,"cd $homedir");}
    if(exists($options->{"-v"})){
        foreach my $key(keys(%{$options->{"-v"}})){
            my $val=$options->{"-v"}->{$key};
            push(@commands,"export $key=$val");
        }
    }
    my $stdoutFile=exists($options->{"-o"})?"$rootdir/".$options->{"-o"}:$stdoutdir;
    mkdir($stdoutFile);
    if(exists($options->{"-N"})){$stdoutFile.="/".$options->{"-N"};}
    mkdir($stdoutFile);
    $stdoutFile.="/$basename.stdout.txt";
    my $stderrFile=exists($options->{"-e"})?"$rootdir/".$options->{"-e"}:$stderrdir;
    mkdir($stderrFile);
    if(exists($options->{"-N"})){$stderrFile.="/".$options->{"-N"};}
    mkdir($stderrFile);
    $stderrFile.="/$basename.stderr.txt";
    $command="$command > $stdoutFile 2> $stderrFile";
    push(@commands,"$command");
    if(exists($options->{"--sleep"})){push(@commands,"sleep ".$options->{"--sleep"});}
    push(@commands,"rm $bashFile");
    push(@commands,"if [ ! -s $stdoutFile ]; then");
    push(@commands,"rm $stdoutFile");
    push(@commands,"fi");
    push(@commands,"if [ ! -s $stderrFile ]; then");
    push(@commands,"rm $stderrFile");
    push(@commands,"fi");
    if(exists($options->{"-N"})){
        push(@commands,"if [ -z \"\$(ls -A $jobDir)\" ]; then");
        push(@commands,"rmdir $jobDir");
        push(@commands,"fi");
    }
    foreach my $line(@commands){print $writer "$line\n";}
    close($writer);
    rename($tmpFile,$bashFile);
    return $bashFile;
}
############################## getDate ##############################
sub getDate{
	my $delim=shift;
	my $time=shift;
	if(!defined($delim)){$delim="";}
	if(!defined($time)||$time eq ""){$time=localtime();}
	else{$time=localtime($time);}
	my $year=$time->year+1900;
	my $month=$time->mon+1;
	if($month<10){$month="0".$month;}
	my $day=$time->mday;
	if($day<10){$day="0".$day;}
	return $year.$delim.$month.$delim.$day;
}
############################## listFilesRecursively ##############################
sub listFilesRecursively{
	my @directories=@_;
	my $filegrep=shift(@directories);
	my $fileungrep=shift(@directories);
	my $recursivesearch=shift(@directories);
	my @inputfiles=();
	foreach my $directory (@directories){
		if(-f $directory){push(@inputfiles,$directory);next;}
		elsif(-l $directory){push(@inputfiles,$directory);next;}
		opendir(DIR,$directory);
		foreach my $file(readdir(DIR)){
			if($file eq "."){next;}
			if($file eq ".."){next;}
			if($file eq ""){next;}
			if($file=~/^\./){next;}
			my $path="$directory/$file";
			if(-d $path){
				if($recursivesearch!=0){push(@inputfiles,listFilesRecursively($filegrep,$fileungrep,$recursivesearch-1,$path));}
				next;
			}
			if(defined($filegrep)&&$file!~/$filegrep/){next;}
			if(defined($fileungrep)&&$file=~/$fileungrep/){next;}
			push(@inputfiles,$path);
		}
		closedir(DIR);
	}
	return sort{$a cmp $b}@inputfiles;
}
############################## openFile ##############################
sub openFile{
	my $path=shift();
	if($path=~/^(.+\@.+)\:(.+)$/){
		if($path=~/\.gz(ip)?$/){return IO::File->new("ssh $1 'gzip -cd $2'|");}
		elsif($path=~/\.bz(ip)?2$/){return IO::File->new("ssh $1 'bzip2 -cd $2'|");}
		elsif($path=~/\.bam$/){return IO::File->new("ssh $1 'samtools view $2'|");}
		else{return IO::File->new("ssh $1 'cat $2'|");}
	}else{
		if($path=~/\.gz(ip)?$/){return IO::File->new("gzip -cd $path|");}
		elsif($path=~/\.bz(ip)?2$/){return IO::File->new("bzip2 -cd $path|");}
		elsif($path=~/\.bam$/){return IO::File->new("samtools view $path|");}
		else{return IO::File->new($path);}
	}
}
############################## printTable ##############################
sub printTable{
	my @out=@_;
	my $return_type=$out[0];
	if(lc($return_type) eq "print"){$return_type=0;shift(@out);}
	elsif(lc($return_type) eq "array"){$return_type=1;shift(@out);}
	elsif(lc($return_type) eq "stderr"){$return_type=2;shift(@out);}
	else{$return_type= 2;}
	printTableSub($return_type,"",@out);
}
sub printTableSub{
	my @out=@_;
	my $return_type=shift(@out);
	my $string=shift(@out);
	my @output=();
	for(@out){
		if(ref( $_ ) eq "ARRAY"){
			my @array=@{$_};
			my $size=scalar(@array);
			if($size==0){
				if($return_type==0){print $string."[]\n";}
				elsif($return_type==1){push(@output,$string."[]");}
				elsif($return_type==2){print STDERR $string."[]\n";}
			}else{
				for(my $i=0;$i<$size;$i++){push(@output,printTableSub($return_type,$string."[$i]=>\t",$array[$i]));}
			}
		} elsif(ref($_)eq"HASH"){
			my %hash=%{$_};
			my @keys=sort{$a cmp $b}keys(%hash);
			my $size=scalar(@keys);
			if($size==0){
				if($return_type==0){print $string."{}\n";}
				elsif($return_type==1){push( @output,$string."{}");}
				elsif($return_type==2){print STDERR $string."{}\n";}
			}else{
				foreach my $key(@keys){push(@output,printTableSub($return_type,$string."{$key}=>\t",$hash{$key}));}
			}
		}elsif($return_type==0){print "$string\"$_\"\n";}
		elsif($return_type==1){push( @output,"$string\"$_\"");}
		elsif($return_type==2){print STDERR "$string\"$_\"\n";}
	}
	return wantarray?@output:$output[0];
}
############################## sortSubs ##############################
sub sortSubs{
	my $path="$prgdir/$program_name";
	my $reader=openFile($path);
	my @headers=();
	my $name;
	my $blocks={};
	my $block=[];
	my $date=getDate("/");
	my @orders=();
	while(<$reader>){
		chomp;s/\r//g;
		if(/^#{30}\s*(\S+)\s*#{30}$/){
			$name=$1;
			if($name!~/^[A-Z]+$/){push(@{$block},$_);last;}
		}elsif(/^my \$program_version=\"\S+\";/){$_="my \$program_version=\"$date\";";}
		push(@headers,$_);
	}
	while(<$reader>){
		chomp;s/\r//g;
		if(/^#{30}\s*(\S+)\s*#{30}$/){
			$blocks->{$name}=$block;
			push(@orders,$name);
			$name=$1;
			$block=[];
		}
		push(@{$block},$_);
	}
	close($reader);
	if(defined($name)){$blocks->{$name}=$block;push(@orders,$name);}
	my ($writer,$file)=tempfile("scriptXXXXXXXXXX",DIR=>".",SUFFIX=>".pl",UNLINK=>1);
	foreach my $line(@headers){print $writer "$line\n";}
	foreach my $key(sort{$a cmp $b}@orders){foreach my $line(@{$blocks->{$key}}){print $writer "$line\n";}}
	close($writer);
    system("mv $file $path");
    system("chmod 755 $path");
}
############################## test ##############################
sub test{
    my ($writer,$tmpFile)=tempfile("bashXXXXXXXXXX",DIR=>"/tmp",SUFFIX=>".sh");
    print $writer "qsub.pl -o log -e log -v first=Hello,second=World -N hello -v number=1 -hold_jid waitID echo \\\$first \\\$second \\\$number\n";
    print $writer "qsub.pl -o log -e log -v first=Hello,second=World -N hello -v number=2 -hold_jid waitID echo \\\$first \\\$second \\\$number\n";
    print $writer "qsub.pl -o log -e log -v first=Hello,second=World -N hello -v number=3 -hold_jid waitID echo \\\$first \\\$second \\\$number\n";
    print $writer "qsub.pl -cwd -N name -sync y -hold_jid hello ls -lt\n";
    system("bash $tmpFile");
}
############################## waitForFreeJobQueue ##############################
sub waitForFreeJobQueue{
	my @files=listFilesRecursively("\.sh",undef,1,$bashdir);
	while(scalar(@files)>=$maxJob){
		sleep($sleepTime);
		@files=listFilesRecursively("\.sh",undef,1,$bashdir);
	}
}
