#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Digest::MD5;
use IPC::System::Simple qw(capture);

# 处理命令行参数
my $output_dir;  # 输出目录
my $srr_accession;  # SRR访问号
GetOptions(
    'output-dir=s' => \$output_dir,
    'srr=s' => \$srr_accession,
);

# 检查参数是否正确提供
unless ($output_dir && $srr_accession) {
    die "Usage: perl download_srr.pl --output-dir <output_directory> --srr <SRR_accession>\n";
}

# 创建输出目录（如果不存在）
mkdir $output_dir unless -d $output_dir;

# 构建SRR下载路径
my $sra_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${srr_accession}/${srr_accession}.sra";

# 下载SRR文件
my $sra_file = "$output_dir/${srr_accession}/${srr_accession}.sra";
system("prefetch", "-O", $output_dir, $srr_accession);

# 检查下载是否成功
unless (-e $sra_file) {
    die "Failed to download SRR file: $sra_file\n";
}

# 计算下载文件的MD5校验和
my $md5_checksum = capture("md5sum", $sra_file);
$md5_checksum =~ s/^\s+|\s+$//g;  # 去除前导和尾随空格
my ($downloaded_md5) = split(' ', $md5_checksum);

# 获取NCBI提供的MD5校验和
my $expected_md5_file = "$sra_file.md5";
my $expected_md5 = "";
if (-e $expected_md5_file) {
    $expected_md5 = capture("cat", $expected_md5_file);
    $expected_md5 =~ s/^\s+|\s+$//g;  # 去除前导和尾随空格
    ($expected_md5) = split(' ', $expected_md5);
} else {
    die "MD5 checksum file not found: $expected_md5_file\n";
}

# 比较两个MD5校验和
if ($downloaded_md5 eq $expected_md5) {
    print "Downloaded file verified successfully.\n";
} else {
    print "Downloaded file verification failed.\n";
}

# 输出SRR文件的绝对路径
print "SRR file downloaded to: $sra_file\n";
