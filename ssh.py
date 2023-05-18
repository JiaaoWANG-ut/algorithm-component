import paramiko
import hashlib
import time
import struct
import hmac
import base64



YourOwnSeries='NRYKB66RYGM75JKJ7T6UANAWCIVRALHP'
refresh: int=30       ###refresh time
token_length: int=6  ###token_length
timestamp = time.time()
counter = int(timestamp) // refresh
msg = struct.pack('>Q', counter)
digest = hmac.new(base64.b32decode(YourOwnSeries), msg, hashlib.sha1).digest()
ob = digest[19]
pos = ob & 15
base = struct.unpack('>I', digest[pos:pos + 4])[0] & 0x7fffffff
token = base % (10**token_length)

#print(str("\n\nToken has copied to your Clipboard.\n\nHave a Nice Day\nJiaao@UT AUSTIN_Henkelman Group\nwangjiaao0720@utexas.edu"))
#print(password)


hostname='login2.frontera.tacc.utexas.edu'
port=22
username='jiaao'
key_filename="id_rsa"
password=str(token)
#password="596263"
#timeout=30

ssh = paramiko.SSHClient()

ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy)
try:
    ssh.connect(hostname=hostname, port=port,
                        username=username, 
                        password=password,
                        key_filename=key_filename,
                        #timeout=timeout,
                        #passphrase=passphrase,
                        compress=True,
                        )
    stdin, stdout, stderr = ssh.exec_command(
        'python3 scripts/nodesview1.1.py;\
            cd /home1/08197/jiaao/Scratch/GCN/trainingset;\
                python3 sbatch.py'
        )
    print(stdout.read().decode())
    ssh.close()
    
except:
    print("Request Busy, try 30s later again")
    ssh.close()




