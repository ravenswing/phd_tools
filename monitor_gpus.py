import pandas as pd
from fabric import Connection

outfile = "/home/rhys/iqtc_today.txt"

# Make the connection to IQTC.
c = Connection(host="portal.qt.ub.edu", user="g19sllabres")

users = ["g19sllabres", "g19torces", "g19fjluque", "g19tviayna"]


def scan_node(n):
    nvidia_out = c.run(f"ssh nodeg{n} -t 'nvidia-smi'", hide="both").stdout.split("\n")
    processes = c.run(
        f"ssh nodeg{n} -t 'ps axfo user:20,pid,start,comm:20'", hide="both"
    ).stdout.split("\n")

    nvidia_out = nvidia_out[
        [i for i, s in enumerate(nvidia_out) if "===" in s][-1] + 1 : -2
    ]
    nvidia_out = [[s.split()[i] for i in (1, 4, 6)] for s in nvidia_out]
    gpus = pd.DataFrame(nvidia_out, columns=["gpu", "pid", "command"])

    processes = [s for s in processes[1:] if s[:4] != "root" and len(s) > 1]
    cut_proc = [s for s in processes if s.split()[1] in gpus.pid.to_list()]
    users = pd.DataFrame(
        [s.replace(" | ", "").replace("\\_ ", "").split() for s in cut_proc],
        columns=["user", "pid", "start", "ps_command"],
    )

    parents = []
    for s in cut_proc:
        if "\\_" in s:
            i = processes.index(s)
            parent = [s]
            while "\\_" in s:
                s = processes[i - 1]
                parent.append(s)
                if "sge_shepherd" in s:
                    break
                else:
                    i -= 1
            parents.append(parent[::-1])

    users["parent"] = parents

    df = pd.merge(gpus, users, on="pid")

    df["node"] = n

    df = df.loc[df.duplicated(subset="gpu", keep=False)]

    return None if df.empty else df


def running_jobs():
    jobs = c.run(
        "export SGE_ROOT=/sge;/sge/bin/lx24-amd64/qstat -u '*' | grep 'iqtc10_g1819_gpu\|fartorgpu'",
        hide="both",
    ).stdout.split("\n")
    jobs = [[s.split()[i] for i in (0, 2, 3, 7)] for s in jobs if len(s) > 1]

    df = pd.DataFrame(jobs, columns=["jobID", "jobname", "user", "node"])
    df.node = df.node.map(lambda x: int(x.split("@nodeg")[-1]))
    return df


def process_df(df, N, jobs):
    msg = []
    for gpu, data in df.groupby("gpu"):
        msg.append(f"\nNodeg{N} - Clash on GPU {gpu}")
        for _, row in data.iterrows():
            msg.append(f"    {row.user:<12} running {row.command}")
            msg.append("    Possible Jobs:")
            possible_jobs = jobs[(jobs.user == row.user) & (jobs.node == N)]
            for _, pj in possible_jobs.iterrows():
                msg.append(f"        {pj.jobID:<8} {pj.jobname}...")

    return msg


def main():
    output = ["13/06/2024"]
    for node in range(11, 17):
        clashes = scan_node(node)
        jobs = running_jobs()

        if clashes is None:
            output.append(f"\nNodeg{node} - No clashes")

        elif not any(user in users for user in clashes.user.to_list()):
            output.append(f"\nNodeg{node} - Clashes but our users not affected.")

        else:
            output.extend(process_df(clashes, node, jobs))

    with open(outfile, "w") as f:
        f.writelines([f"{line}\n" for line in output])


if __name__ == "__main__":
    main()
