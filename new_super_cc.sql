select UNIX_TIMESTAMP(record_ts,'yyyy-MM-dd-hh.mm') as epoch, cidr, 
COALESCE(sum(a.group_map['Automotive']),0) as Automotive,
COALESCE(sum(a.group_map['Business Services']),0) as `Business Services`,
COALESCE(sum(a.group_map['High Technology']),0) as `High Technology`,
COALESCE(sum(a.group_map['Media & Entertainment']),0) as `Media & Entertainment`,
COALESCE(sum(a.group_map['Software as a Service']),0) as `Software as a Service`,
COALESCE(sum(a.group_map['Hotel & Travel']),0) as `Hotel & Travel`,
COALESCE(sum(a.group_map['Foundation-Not for Profit']),0) as `Foundation-Not for Profit`,
COALESCE(sum(a.group_map['Manufacturing']),0) as Manufacturing,
COALESCE(sum(a.group_map['Financial Services']),0) as `Financial Services`,
COALESCE(sum(a.group_map['Consumer Goods']),0) as `Consumer Goods`,
COALESCE(sum(a.group_map['Consumer Services']),0) as `Consumer Services`,
COALESCE(sum(a.group_map['Energy & Utilities']),0) as `Energy & Utilities`,
COALESCE(sum(a.group_map['Pharma/Health Care']),0) as `Pharma/Health Care`,
COALESCE(sum(a.group_map['Public Sector']),0) as `Public Sector`,
COALESCE(sum(a.group_map['Energy & Utilities']),0) as `Energy & Utilities`,
COALESCE(sum(a.group_map['Real Estate']),0) as `Real Estate`,
COALESCE(sum(a.group_map['Education']),0) as Education,
COALESCE(sum(a.group_map['Gaming']),0) as Gaming,
COALESCE(sum(a.group_map['Gambling']),0) as Gambling
from
( select record_ts, cidr, 
map(super_cat,suc_hits) as group_map
from tosi7_072017_super_cc
where cidr in ('104.218.24.0/24', '12.10.219.0/24', '12.163.3.0/24', '12.165.188.0/24', '12.216.81.0/24', '12.96.67.0/24', '129.83.31.0/24', '13.124.166.0/24', '13.124.179.0/24', '131.242.9.0/24', '132.245.35.0/24', '132.245.65.0/24', '132.245.76.0/24', '134.132.40.0/24', '134.170.68.0/24', '13.55.220.0/24', '136.147.62.0/24', '136.200.53.0/24', '136.32.184.0/24', '136.57.5.0/24', '136.62.170.0/24', '136.62.94.0/24', '138.68.130.0/24', '138.68.137.0/24', '138.68.43.0/24', '139.60.64.0/24', '139.99.130.0/24', '143.115.159.0/24', '143.121.239.0/24', '144.212.3.0/24', '145.8.179.0/24', '146.235.66.0/24', '146.235.89.0/24', '146.241.163.0/24', '146.241.208.0/24', '146.241.33.0/24', '146.241.48.0/24', '15.203.233.0/24', '152.51.48.0/24', '155.109.35.0/24', '155.186.118.0/24', '155.186.82.0/24', '155.192.33.0/24', '155.254.17.0/24', '157.7.223.0/24', '159.203.218.0/24', '159.214.124.0/24', '161.141.6.0/24', '161.150.176.0/24', '161.49.249.0/24', '162.119.11.0/24', '162.210.14.0/24', '162.211.75.0/24', '162.254.24.0/24', '163.166.8.0/24', '164.85.24.0/24', '165.20.108.0/24', '165.204.77.0/24', '167.206.147.0/24', '169.198.254.0/24', '170.109.16.0/24', '170.128.128.0/24', '170.131.131.0/24', '170.194.16.0/24', '170.22.76.0/24', '170.236.180.0/24', '170.249.184.0/24', '170.40.193.0/24', '170.61.236.0/24', '17.132.107.0/24', '17.142.142.0/24', '17.224.126.0/24', '17.253.39.0/24', '17.255.235.0/24', '173.161.78.0/24', '173.164.81.0/24', '173.165.201.0/24', '173.224.161.0/24', '173.227.210.0/24', '173.244.130.0/24', '176.10.104.0/24', '182.16.234.0/24', '182.171.248.0/24', '183.82.0.0/24', '185.31.7.0/24', '185.63.117.0/24', '188.68.0.0/24', '190.4.3.0/24', '192.136.15.0/24', '192.158.48.0/24', '192.162.61.0/24', '192.30.38.0/24', '193.127.224.0/24', '193.99.214.0/24', '194.120.84.0/24', '194.126.214.0/24', '194.60.38.0/24', '195.110.84.0/24', '195.128.10.0/24', '195.19.82.0/24', '195.7.17.0/24', '198.169.6.0/24', '198.175.53.0/24', '198.176.84.0/24', '198.177.1.0/24', '198.186.162.0/24', '198.199.206.0/24', '198.22.21.0/24', '198.233.245.0/24', '198.28.69.0/24', '198.36.40.0/24', '199.101.6.0/24', '199.244.219.0/24', '199.248.185.0/24', '199.250.80.0/24', '199.87.0.0/24', '199.96.153.0/24', '200.185.183.0/24', '20.138.2.0/24', '202.161.57.0/24', '202.214.112.0/24', '202.232.184.0/24', '202.232.192.0/24', '202.238.159.0/24', '202.4.248.0/24', '203.110.235.0/24', '204.115.199.0/24', '204.145.112.0/24', '204.152.239.0/24', '204.68.32.0/24', '204.76.134.0/24', '204.87.189.0/24', '205.144.216.0/24', '206.200.255.0/24', '206.223.21.0/24', '207.144.152.0/24', '207.156.9.0/24', '208.79.246.0/24', '209.102.210.0/24', '212.31.118.0/24', '212.33.148.0/24', '212.56.120.0/24', '212.65.8.0/24', '213.68.126.0/24', '214.28.226.0/24', '216.10.224.0/24', '216.115.13.0/24', '216.136.78.0/24', '216.52.85.0/24', '217.110.243.0/24', '217.23.230.0/24', '217.89.158.0/24', '24.97.34.0/24', '34.194.79.0/24', '34.207.222.0/24', '34.207.53.0/24', '34.227.113.0/24', '34.228.170.0/24', '34.228.17.0/24', '34.239.249.0/24', '38.117.92.0/24', '38.88.200.0/24', '40.100.150.0/24', '40.101.64.0/24', '40.101.72.0/24', '40.103.47.0/24', '40.128.64.0/24', '40.96.21.0/24', '40.97.115.0/24', '40.97.132.0/24', '40.97.159.0/24', '40.97.165.0/24', '40.97.171.0/24', '41.203.114.0/24', '45.64.99.0/24', '46.101.45.0/24', '46.165.168.0/24', '46.34.80.0/24', '47.21.131.0/24', '50.112.119.0/24', '50.19.124.0/24', '50.192.75.0/24', '50.203.182.0/24', '50.243.86.0/24', '50.73.49.0/24', '50.76.77.0/24', '5.2.198.0/24', '52.204.211.0/24', '52.205.102.0/24', '52.208.107.0/24', '52.211.105.0/24', '52.221.181.0/24', '52.54.209.0/24', '52.54.92.0/24', '52.55.47.0/24', '58.227.74.0/24', '58.26.32.0/24', '61.206.119.0/24', '63.157.110.0/24', '63.247.60.0/24', '63.78.207.0/24', '64.109.211.0/24', '64.117.234.0/24', '64.132.0.0/24', '64.206.241.0/24', '64.53.145.0/24', '65.116.116.0/24', '65.122.111.0/24', '65.210.80.0/24', '65.90.11.0/24', '66.35.55.0/24', '66.76.179.0/24', '67.132.198.0/24', '67.20.180.0/24', '67.223.1.0/24', '67.59.68.0/24', '68.168.104.0/24', '69.167.206.0/24', '69.17.145.0/24', '69.196.152.0/24', '69.4.98.0/24', '69.57.50.0/24', '70.104.26.0/24', '72.22.31.0/24', '72.36.1.0/24', '72.91.227.0/24', '74.117.207.0/24', '75.126.38.0/24', '80.169.129.0/24', '80.169.158.0/24', '80.245.147.0/24', '80.85.85.0/24', '81.110.143.0/24', '81.173.209.0/24', '81.20.89.0/24', '82.147.18.0/24', '87.118.116.0/24', '91.102.14.0/24', '91.194.221.0/24', '96.47.146.0/24', '96.63.178.0/24', '97.64.197.0/24', '98.190.206.0/24') a 
group by a.record_ts, a.cidr
order by cidr asc, epoch asc
limit 225000;

