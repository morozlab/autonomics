import argparse
import threading
from SimpleXMLRPCServer import SimpleXMLRPCServer
from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler
from zeroclick import settings

class RequestHandler(SimpleXMLRPCRequestHandler):
    rpc_paths = ('/RPC2')


max_resources = {}
resources = {}


def add(d, location):
    if(not location in max_resources):
        max_resources[location] = {}
        resources[location] = {}

    this_res, this_max = _get_res_dicts(location)
    for res, amt in d.items():
        if(res in this_max):
            this_max[res] += amt
            this_res[res] += amt
                
        else:
            this_max[res] = amt
            this_res[res] = amt

    return this_max


def current_resources():
    return (resources, max_resources)


def _get_res_dicts(location):
    return (resources[location], max_resources[location])

            
def relinquish(d, location):
    this_res, this_max = _get_res_dicts(location)
    for res, amt in d.items():
        if(res in this_res):
            if(this_res[res] + amt > this_max[res]):
                this_res[res] = this_max[res]
            else:
                this_res[res] += amt

    return this_res


def remove(d, location):
    this_res, this_max = _get_res_dicts(location)
    for res, amt in d.items():
        if(res in this_max):
            this_max[res] -= amt
            if(this_max[res] < 0):
                this_max[res] = 0

    return this_max


def request(d, location):
    this_res = _get_res_dicts(location)[0]
    for res, amt in d.items():
        if(res in this_res):
            if(amt > this_res[res]):
                return False

    for res, amt in d.items():
        if(res in this_res):
            this_res[res] -= amt

    return True


def reset_to_max(location):
    this_res, max_res = _get_res_dicts(location)
    this_res = {key: value for key, value in max_res.items()}
    return True

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--resources', dest='resources', nargs='+',
                        default=None,
                        help='''A space-separated list
                                of starting resources for
                                this resource server. Ex:
                                location:cpu:50 location:pfam:1
                                location:blast_nr:1
                                ''')

    args = parser.parse_args()

    if(not args.resources is None):
        for resource in args.resources:
            els = resource.split(":")
            location = els[0]
            if(not location in max_resources):
                max_resources[location] = {}
                resources[location] = {}
            resource = els[1]
            amt = int(els[2])
            max_resources[location][resource] = amt 
            resources[location][resource] = amt
            
    s = SimpleXMLRPCServer(('localhost',
                            settings.RESOURCE_SERVER_PORT),
                            requestHandler=RequestHandler)

    s.register_introspection_functions()
    s.register_function(add)
    s.register_function(current_resources)
    s.register_function(remove)
    s.register_function(relinquish)
    s.register_function(request)
    s.serve_forever()
    
