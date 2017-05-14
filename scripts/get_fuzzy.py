#!/Users/ponet/tools/anaconda3/envs/PhD/bin/python
if __name__ == "__main__":
    import sys
    from DFTools import SSHSession

    sysarg = sys.argv
    ext_dir = sysarg[1]
    loc_dir = sysarg[2]
    fuzzy = sysarg[3]

    usr = "ponet"
    host = "10.255.9.115"

    ssh=SSHSession(host,usr,key_file=open("/Users/ponet/.ssh/id_rsa"))
    ssh.get_fuzzy(ext_dir,loc_dir,fuzzy)
    print("done")
    
