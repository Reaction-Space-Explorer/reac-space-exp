import os

def call_neo4j(network):
    os.execute(f"python neo4j_loader_and_queries/main.py {network}")

# choose network to generate
argument = "glucose"
if argument == "glucose":
    # glucose network
    os.execute("mod -f main/main.py")
    call_neo4j(network = "glucose")
elif argument == "radicals":
    # radicals network
    os.execute("mod -f radicals/all7/all7.py")

