import pickle

def save_obj(obj, name, directory):
    with open(directory + '/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name, directory):
    with open(directory + '/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
