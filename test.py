# Set environment workspace in current directory
path_script = os.path.dirname(os.path.abspath(__file__))
env.workspace = os.path.join(path_script, 'demo_data', 'benchmark_data')

#env.workspace = r"E:\Mix\WYZ_AcademicWriting\paper_pycor\arcPycor\demo_data\benchmark_data"

# Folder for outputs
dirOutputs = os.path.join(env.workspace, 'outputs')
if os.path.exists(dirOutputs):
    shutil.rmtree(dirOutputs)
os.makedirs(dirOutputs)