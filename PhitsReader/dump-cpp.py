import ctypes
import json
import os
import shutil

# JSONファイルの読み込み
with open("G:/tagawa/phits/trace/input.json", "r") as f:
    para = json.load(f)

# ファイルパスを文字列フォーマットで指定
DatPath = f"G:/tagawa/phits/trace/output_8/{para['output']}/dumpall.dat".encode('utf-8')
OutputDirectory = f"G:/tagawa/phits/trace/output_8/{para['output']}"

# DLLの読み込み
PhitsReader = ctypes.CDLL('./PhitsReader.dll')

# 関数の引数と戻り値の型を指定
PhitsReader.MakeOutput.argtypes = [ctypes.c_char_p]
PhitsReader.MakeOutput.restype = None

# DLL関数の呼び出し
PhitsReader.MakeOutput(DatPath)

CurrentDirectory = os.path.dirname(os.path.abspath(__file__))
OutputFile = 'output.json'

# 移動先ディレクトリが存在しない場合、作成
if not os.path.exists(OutputDirectory):
    os.makedirs(OutputDirectory)

# ファイルのパスを構築
OutputCurrentPath = os.path.join(CurrentDirectory, OutputFile)
OutputDestinationPath = os.path.join(OutputDirectory, OutputFile)

# ファイルが存在するかを確認
if os.path.isfile(OutputCurrentPath):
    try:
        shutil.move(OutputCurrentPath, OutputDestinationPath)
        print(f'Moved {OutputFile} to {OutputDirectory}')
    except Exception as e:
        print(f'Error moving {OutputFile}: {e}')
else:
    print(f'{OutputFile} does not exist in the current directory')
